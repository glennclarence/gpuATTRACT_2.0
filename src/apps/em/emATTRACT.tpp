/*
 * mcATTRACT.tpp
 *
 *  Created on: Jul 11, 2017
 *      Author: uwe
 */

#ifndef SRC_APPS_EM_EMATTRACT_TPP_
#define SRC_APPS_EM_EMATTRACT_TPP_

#include <iostream>
#include <fstream>
#include "emATTRACT.h"
#include "Configurator_6D.h"
#include "Request.h"
#include "Server.h"
#include "RequestHandler.h"
#include <chrono>
#include "SolverBase.h"


namespace as {

template<typename GenericTypes>
emATTRACT<GenericTypes>::emATTRACT() : _config(new configurator_t()) {}

template<typename GenericTypes>
void emATTRACT<GenericTypes>::init(CmdArgs const& args) {
	_config->init(args);
	trackGradients = args.trackGradients;
	trackStates = args.trackStates;
	minModesOnly = args.minmodesonly;
	track_file = args.trackFile;
	mode_thresh = args.mode_thresh;
}

template<typename GenericTypes>
void emATTRACT<GenericTypes>::finalize() {
	_config->finalize();
}

template<typename GenericTypes>
void emATTRACT<GenericTypes>::run() {

	auto& dofs = _config->dofs();
	auto server = _config->server();
	auto& common = _config->common();

	RequestHandler<GenericTypes> requestHandler = RequestHandler<GenericTypes>::newBuilder()
			.withServer(server)
			.withCommon(common)
			.withDofs(dofs)
			.withSolverName("VA13")
			.withSolverSettings(trackStates,trackGradients, minModesOnly, mode_thresh)
			.build();
	//std::cout << " track grads " << trackGradients << " trackStates " << trackStates << " min modes only " << minModesOnly << std::endl;
	auto start = std::chrono::system_clock::now();
	requestHandler.run();
	auto end = std::chrono::system_clock::now();
// THis is just for testing purposes
	auto results_grad = requestHandler.getResultEnGrads();
	auto results_dof = requestHandler.getResultStates();
	int count = 0;
	//for (auto const res : results_dof) {


//	std::cout << "#pivot 1 "<< " "<< 8.95646 << " "<<56.9144 << " "<<92.5237 << std::endl;
//	std::cout << "#pivot 2 "<< " "<<8.91543 << " "<<34.5702 << " "<<114.345 << std::endl;

	std::cout << "#pivot 1 "<< " "<< common.pivotRec.x << " "<<common.pivotRec.y << " "<<common.pivotRec.z << std::endl;
	std::cout << "#pivot 2 "<< " "<<common.pivotLig.x << " "<<common.pivotLig.y << " "<<common.pivotLig.z << std::endl;
	std::cout << "#centered receptor: false "<< std::endl;
	std::cout << "#centered ligand: false "<< std::endl;
	for (int i = 0; i < results_dof.size(); i++){
		std::cout << "#"<< i+1<< std::endl;
		std::cout <<"## Energy: " << results_grad[i].get_Energy() << std::endl;
		auto  pos = results_dof[i].get_pos();
				pos.x = pos.x + common.pivotRec.x -  common.pivotLig.x;
				pos.y = pos.y + common.pivotRec.y -  common.pivotLig.y;
				pos.z = pos.z + common.pivotRec.z -  common.pivotLig.z;
				results_dof[i].set_pos(pos.x,pos.y, pos.z);
				std::cout << results_dof[i] << std::endl;


	}

	/*std::vector<std::shared_ptr<std::vector<std::vector<float>>>> ts = requestHandler.getResultStateTracker();
	for (int i = 0; i < ts.size(); i++){
	for (int k = 0; k < ts[i]->size(); k++){
		auto ref = *ts[i];
		for (int l = 0; l < ref[k].size(); l++){
		std::cout << ref[k][l] <<" ";
	}
	std::cout<< std::endl;
	}
		std::cout<< std::endl;
	}


	std::vector<std::shared_ptr<std::vector<std::vector<float>>>> ts2 = requestHandler.getResultEnGradTracker();
		for (int i = 0; i < ts2.size(); i++){
		for (int k = 0; k < ts2[i]->size(); k++){
			auto ref = *ts2[i];
			for (int l = 0; l < ref[k].size(); l++){
			std::cout << ref[k][l] <<" ";
		}
		std::cout<< std::endl;
		}
			std::cout<< std::endl;
		}*/

		if( trackStates || trackGradients ){

			std::vector<std::shared_ptr<std::vector<std::vector<float>>>> gradients = requestHandler.getResultEnGradTracker();
			std::vector<std::shared_ptr<std::vector<std::vector<float>>>> states = requestHandler.getResultStateTracker();

			auto track = requestHandler.getTrack();
			std::fstream fs;
			fs.open (track_file, std::fstream::in | std::fstream::out | std::fstream::trunc );
			fs <<"dofIdx,step";
			if(trackGradients){
				fs << ",energy";
				if(minModesOnly == 0){
						fs << ",g_rot_x,g_rot_y,g_rot_z,g_trans_x,g_trans_y,g_trans_z";
				}
				for( int i = 0; i< Common_Modes::numModesRec; ++i){
					fs <<",g_mode_rec_"<< i+1;
				}
				for( int i = 0; i< Common_Modes::numModesLig; ++i){
					fs <<",g_mode_lig_"<< i+1;
				}



			}
			if(trackStates){
				if(minModesOnly == 0){
				fs << ",s_rot_x,s_rot_y,s_rot_z,s_trans_x,s_trans_y,s_trans_z";
				}
				for( int i = 0; i< Common_Modes::numModesRec; ++i){
					fs <<",s_mode_rec_"<< i+1;
				}
				for( int i = 0; i< Common_Modes::numModesLig; ++i){
					fs <<",s_mode_lig_"<< i+1;
				}
			}
			fs <<std::endl;
			for (int i = 0; i < track.size(); i++){
				for (int k = 0; k < track[i]->size(); k++){
					auto ref = *track[i];
					fs <<i+1<<","<<k;
					if(trackGradients){
					fs<<"," << ref[k].energy;
					for ( auto g : ref[k].grads){
						fs<<"," << g ;
					}}
					if(trackStates){
					for ( auto s : ref[k].states){
						fs<<"," << s;
					}
					}
					fs<< std::endl;
				}
				}


			fs.close();
		}

		//std::cout << "time"<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << std::endl;
}

} //namespace



#endif /* SRC_APPS_EM_EMATTRACT_TPP_ */
