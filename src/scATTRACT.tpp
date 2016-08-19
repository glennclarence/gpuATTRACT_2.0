/*
 * scATTRACT.tpp
 *
 *  Created on: Aug 17, 2016
 *      Author: uwe
 */

#ifndef SRC_SCATTRACT_TPP_
#define SRC_SCATTRACT_TPP_

#include "scATTRACT.h"
#include "Configurator_6D.h"
#include "Request.h"
#include "Server.h"

namespace as {

template<typename SERVICE>
scATTRACT<SERVICE>::scATTRACT() : _config(new configurator_t()) {}

template<typename SERVICE>
void scATTRACT<SERVICE>::init(CmdArgs const& args) {
	_config->init(args);
}

template<typename SERVICE>
void scATTRACT<SERVICE>::finalize() {
	_config->finalize();
}

template<typename SERVICE>
void scATTRACT<SERVICE>::run() {

	auto& dofs = _config->dofs();
	auto& server = _config->server(); //this->server();
	auto& common = _config->common();
	size_t numDofs = dofs.size();
	auto results = std::vector<result_t>(dofs.size());
	Request<dof_t, common_t> request(dofs.data(), numDofs, common);
	server.submit(request);


}

}  // namespace as



#endif /* SRC_SCATTRACT_TPP_ */
