/*
 * GPU_6D_EnergyService.tpp
 *
 *  Created on: Sep 8, 2016
 *      Author: uwe
 */
/*
 * GPU_6D_EnergyService.tpp
 *
 *  Created on: Sep 8, 2016
 *      Author: uwe
 */

#ifndef SRC_SERVICE_GPUENERGYSERVICE6DMODES_TPP_
#define SRC_SERVICE_GPUENERGYSERVICE6DMODES_TPP_

#ifdef CUDA

#include <nvToolsExt.h>
#include "GPUEnergyService6DModes.h"

#include <cassert>
#include "WorkerBuffer.h"
#include "DataManager.h"
#include "Allocator.h"
#include "RingArray.h"
#include "DataItem.h"
#include "DeviceItem.h"
#include "WorkItem.h"
#include "DeviceProtein.h"
#include "DeviceGridUnion.h"
#include "DeviceParamTable.h"
#include "SimParam.h"

#include "transform_modes.h"
#include "interpolation.h"
#include "neighborlist_modes.h"
#include "reduction_modes.h"
#include "macros.h"

#include <iomanip>
#include <limits>
#include <mutex>


namespace as {

template<typename REAL>
GPUEnergyService6DModes<REAL>::GPUEnergyService6DModes(std::shared_ptr<DataManager> dataMng,
		std::vector<int> const& deviceIds) :
	GPUEnergyService<Types_6D_Modes<REAL>>::GPUEnergyService(dataMng), _workerId(0), _deviceIds(deviceIds), _threadsPerDevice(1)
{}

template<typename REAL>
GPUEnergyService6DModes<REAL>::GPUEnergyService6DModes(std::shared_ptr<DataManager> dataMng,
		std::vector<int> const& deviceIds, uint threadsPerDevice) :
	GPUEnergyService<Types_6D_Modes<REAL>>::GPUEnergyService(dataMng), _workerId(0), _deviceIds(deviceIds), _threadsPerDevice(threadsPerDevice)
{}

template<typename REAL>
struct GPUEnergyService6DModes<REAL>::StageResource {
private:
	using workItem_t = typename GPUEnergyService6DModes<REAL>::workItem_t;
public:
	d_GridUnion<REAL> gridRec;
	d_GridUnion<REAL> gridLig;
	d_Protein<REAL>* rec;
	d_Protein<REAL>* lig;
	d_ParamTable<REAL>* table;
	SimParam<REAL>* simParam;
	workItem_t* item;
};

template<typename REAL>
auto GPUEnergyService6DModes<REAL>::createStageResource(workItem_t* item, unsigned const& deviceId) -> StageResource {
	const auto common = item->common();
	auto gridRec = std::dynamic_pointer_cast<DeviceGridUnion<REAL>>(this->_dataMng->get(common->gridIdRec, deviceId)).get(); // _dataMng->get() returns shared_ptr<DataItem>
	assert(gridRec != nullptr);

	auto gridLig = std::dynamic_pointer_cast<DeviceGridUnion<REAL>>(this->_dataMng->get(common->gridIdLig, deviceId)).get(); // _dataMng->get() returns shared_ptr<DataItem>
	assert(gridLig != nullptr);

	auto rec = std::dynamic_pointer_cast<DeviceProtein<REAL>>(this->_dataMng->get(common->recId, deviceId)).get();
	assert(rec != nullptr);

	auto lig = std::dynamic_pointer_cast<DeviceProtein<REAL>>(this->_dataMng->get(common->ligId, deviceId)).get();
	assert(lig != nullptr);

	auto table = std::dynamic_pointer_cast<DeviceParamTable<REAL>>(this->_dataMng->get(common->tableId, deviceId)).get();
	assert(table != nullptr);

	auto simParam = std::dynamic_pointer_cast<SimParam<REAL>>(this->_dataMng->get(common->paramsId)).get();
	assert(simParam != nullptr);

	StageResource stageResource;
	stageResource.gridRec 	= gridRec->getDesc();
	stageResource.gridLig 	= gridLig->getDesc();
	stageResource.lig 		= &lig->desc;
	stageResource.rec 		= &rec->desc;
	stageResource.table 	= &table->desc;
	stageResource.simParam = simParam;
	stageResource.item 	= item;

	return stageResource;
}

template<typename REAL>
class GPUEnergyService6DModes<REAL>::Private {
	using dof_t = typename GPUEnergyService6DModes<REAL>::input_t;
	using workItem_t = typename GPUEnergyService6DModes<REAL>::workItem_t;

public:

	Private(int deviceId)  {
        int i1d;
        cudaSetDevice(deviceId);
        cudaGetDevice(&i1d);
        CUDA_CHECK(cudaStreamCreate(&_stream));
	}
	/**
	 * Allocate new Buffers with size. Old buffers are automatically deallocated;
	 */
	void allocateBuffer( size_t const& numDOFs, size_t const& numAtomsRec, size_t const& numAtomsLig, size_t const& dofSizeRec, size_t const& dofSizeLig ) {
		const size_t atomBufferSizeRec = numDOFs*numAtomsRec;
		const size_t atomBufferSizeLig = numDOFs*numAtomsLig;

		d_defoRec = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(3,atomBufferSizeRec));
		d_defoLig = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(3,atomBufferSizeLig));
		d_potRec = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(4,atomBufferSizeRec));
		d_potLig = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(4,atomBufferSizeLig));
		d_resRec = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(1, dofSizeRec*numDOFs));
		d_resLig = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(1, dofSizeLig*numDOFs));
		h_resRec = std::move(WorkerBuffer<REAL, HostPinnedAllocator<REAL>>(1, dofSizeRec*numDOFs));
		h_resLig = std::move(WorkerBuffer<REAL, HostPinnedAllocator<REAL>>(1, dofSizeLig*numDOFs));
		d_dof    = std::move(WorkerBuffer<dof_t,DeviceAllocator<dof_t>>(1,numDOFs));
		d_trafoRec = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(3,atomBufferSizeRec));
		d_trafoLig = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(3,atomBufferSizeLig));
	}

	size_t bufferSize( ) const {
		return d_trafoLig.bufferSize(  );
	}

	void addItemAndLigandSize(StageResource const& resc) {
		_resources = resc;
	}

	void resizeBuffersIfRequired(size_t const& numDOFs, size_t const& numAtomsRec, size_t const& numAtomsLig, size_t const& dofSizeRec, size_t const& dofSizeLig) {
		if (numDOFs*numAtomsLig > d_trafoLig.bufferSize() || numDOFs*numAtomsRec > d_trafoRec.bufferSize()) {
			allocateBuffer(numDOFs, numAtomsRec, numAtomsLig, dofSizeRec, dofSizeLig );
		}
	}

	size_t getSharedMemSize(int const& id){
		cudaDeviceProp deviceProp;
		cudaVerify(cudaGetDeviceProperties(&deviceProp, id));
		return deviceProp.sharedMemPerBlock;
	}

	size_t getMaxBockSize( int const& id, unsigned const dofSize ){
		size_t sharedMem = getSharedMemSize( id );
		size_t pow2 = 2;
		while (pow2* dofSize *sizeof(REAL) < sharedMem) {
			pow2 *= 2;
		}
		size_t blockSize = pow2 / 2;
		return blockSize;
	}

	void configureDevice() {

		if ( blockSizeReduceRec == 0 ) {
			int id;
			cudaVerify(cudaGetDevice(&id));
			blockSizeReduceRec = getMaxBockSize( id, dofSizeRec );
		}
		if ( blockSizeReduceLig == 0 ) {
			int id;
			cudaVerify(cudaGetDevice(&id));
			blockSizeReduceLig = getMaxBockSize( id, dofSizeLig );
		}

	}

	void process_item( ) {
		auto* const it = _resources.item;

		cudaVerify(cudaMemcpyAsync(d_dof.get(0), it->inputBuffer(), it->size()*sizeof(dof_t), cudaMemcpyHostToDevice, _stream));

		const unsigned numAtomsReceptor = it->size()*_resources.rec->numAtoms;
		const unsigned numAtomsLigand = it->size()*_resources.lig->numAtoms;

		assert( numAtomsReceptor <= d_trafoRec.bufferSize() );
		assert( numAtomsLigand <= d_trafoLig.bufferSize() );

		size_t gridSizeRec = ( numAtomsReceptor + BLSZ_TRAFO - 1) / BLSZ_TRAFO;
		size_t gridSizeLig = ( numAtomsLigand + BLSZ_TRAFO - 1) / BLSZ_TRAFO;

		it->setProcessed();
	}



	WorkerBuffer<dof_t, DeviceAllocator<dof_t>> d_dof;
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_defoRec;
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_defoLig;
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_trafoRec;
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_trafoLig;
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_potRec;
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_potLig;
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_resRec;
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_resLig;
	WorkerBuffer<REAL, HostPinnedAllocator<REAL>> h_resRec;
	WorkerBuffer<REAL, HostPinnedAllocator<REAL>> h_resLig;


	static constexpr size_t BLSZ_TRAFO = 128;
	static constexpr size_t BLSZ_INTRPL = 128;
	size_t blockSizeReduceRec = 0;
	size_t blockSizeReduceLig = 0;
	size_t dofSizeRec = 13 + Common_Modes::numModesRec;
	size_t dofSizeLig = 13 + Common_Modes::numModesLig;
	StageResource _resources;
	cudaStream_t _stream; /** cuda streams */
};

template<typename REAL>
auto GPUEnergyService6DModes<REAL>::createDistributor() -> distributor_t {
	distributor_t fncObj = [this] (common_t const* common, size_t numWorkers) {
		std::vector<id_t> ids = {common->gridIdRec, common->gridIdLig, common->ligId, common->recId, common->tableId};
		auto id = this->_dataMng->getCommonDeviceIds(ids);
		std::vector<as::workerId_t> vec(numWorkers);
		std::iota(vec.begin(), vec.end(), 0);
		return vec;
	};
	return fncObj;
}

template<typename REAL>
auto GPUEnergyService6DModes<REAL>::createItemProcessor() -> itemProcessor_t {


	deviceId_t deviceId = _deviceIds[_workerId/_threadsPerDevice];
	std::shared_ptr<Private> p = std::make_shared<Private>(deviceId);
	 _workerId++;
	itemProcessor_t fncObj = [this, deviceId, p] (workItem_t* item) -> bool {

		/* Set the device to work with */
		cudaVerify(cudaSetDevice(deviceId));
		p->configureDevice();
		if (item != nullptr) {

			auto dI = createStageResource(item, deviceId);

			const auto itemSize = item->size();
			assert(itemSize > 0);

			const unsigned& numAtomsRec = dI.rec->numAtoms;
			unsigned const& numAtomsLig = dI.lig->numAtoms;
			unsigned const& dofSizeRec = 13 + Common_Modes::numModesRec;
			unsigned const& dofSizeLig = 13 + Common_Modes::numModesLig;

			p->addItemAndLigandSize(dI);
			p->resizeBuffersIfRequired( itemSize, numAtomsRec, numAtomsLig, dofSizeRec, dofSizeLig);
		} else {
			return false;
		}
		p->process_item();

		return true;

	};

	return fncObj;
}

} // namespace as

#endif

#endif /* SRC_SERVICE_GPUENERGYSERVICE6D_TPP_ */
