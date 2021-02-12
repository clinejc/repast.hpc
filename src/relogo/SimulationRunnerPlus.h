/*
 *   Repast for High Performance Computing (Repast HPC)
 *
 *   Copyright (c) 2010 Argonne National Laboratory
 *   All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with
 *   or without modification, are permitted provided that the following
 *   conditions are met:
 *
 *     Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *
 *     Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *
 *     Neither the name of the Argonne National Laboratory nor the names of its
 *     contributors may be used to endorse or promote products derived from
 *     this software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *   PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE TRUSTEES OR
 *   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 *   EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 *  SimulationRunnerPlus.h
 *
 *  Created on: Sep 30, 2014
 *      Author: Jon Cline <jon.c.cline@gmail.com>
 */

#ifndef SIMULATIONRUNNERPLUS_HPP_
#define SIMULATIONRUNNERPLUS_HPP_

#include <vector>

#include <boost/mpi/communicator.hpp>

#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"

#include "boost/scoped_ptr.hpp"

#include "repast_hpc/Properties.h"
#include "repast_hpc/RepastProcess.h"

#include "RelogoLink.h"
#include "WorldCreator.h"
#include "WorldDefinition.h"

namespace repast {

  namespace relogo {

    using namespace boost;
    using namespace boost::algorithm;

    /**
     * Runs a Relogo simulation.
     *
     * @tparam ObserverType the type of Observer to create. This type must
     * extend relogo::Observer.
     * @tparam PatchType the type of Patches to create. This must extend
     * relogo::Patch.
     */
    template<typename ObserverType, typename PatchType>
    class SimulationRunnerPlus {
    protected:
      boost::mpi::communicator* comm;
      std::unique_ptr<Properties> mCp_props;
      std::unique_ptr<ObserverType> mCp_obs;

    public:
      /**
       * Creates a SimulationRunner.
       *
       * @param world handle to mpi communication
       */
      SimulationRunnerPlus(boost::mpi::communicator* communicator): comm(communicator) {
	mCp_props.reset(new Properties());
      }

      /**
       * Creates the simulation using the properties defined in props.
       * The properties file must have the following properties defined:
       *
       * <ul>
       * <li>min.x the minimum integer x coordinate of the world
       * <li>min.y the minimum integer y coordinate of the world
       * <li>max.x the maximum integer x coordinate of the world
       * <li>max.h the maximum integer y coordinate of the world
       * <li>grid.buffer the size of the grid and space buffers
       * <li>proc.per.x the number of processes to assign to the world's x dimension. proc.per.x
       * multiplied by proc.per.y must equal the number processes that the simulation will run on
       * <li>proc.per.y the number of processes to assign to the world's y dimension. proc.per.x
       * multiplied by proc.per.y must equal the number processes that the simulation will run on
       * <li>stop.at the tick at which to stop the simulation
       *
       * This will create an Observer of the specified type and populate the world
       * with Patches of the specified type. It will then call setup(props) on that Observer
       * implementation and start the simulation schedule which will call the Observer's
       * go method each tick.
       *
       * @param props a properties file containing the properties mentioned above
       */
      void setup(Properties& props);

      /**
       * Runs the simulation.
       */
      void run();

      /**
       * @returns reference to properties (const)
       */
      const repast::Properties& getProperties() {
	return *mCp_props;
      }

      /**
       * @returns handle to Observer
       */
      ObserverType* getObserver() {
	return mCp_obs.get();
      }
    };

    template<typename ObserverType, typename PatchType>
    void SimulationRunnerPlus<ObserverType,PatchType>::setup(Properties& props) {
      *mCp_props = props;

      int minX = strToInt(mCp_props->getProperty("min.x"));
      int minY = strToInt(mCp_props->getProperty("min.y"));
      int maxX = strToInt(mCp_props->getProperty("max.x"));
      int maxY = strToInt(mCp_props->getProperty("max.y"));
      int buffer = strToInt(mCp_props->getProperty("grid.buffer"));

      bool worldIsWrapped = true;
      if(mCp_props->contains("non.toroidal")) worldIsWrapped = false;

      RelogoLinkContentManager rlcm;

      WorldDefinition def(minX, minY, maxX, maxY, worldIsWrapped, buffer);

      // Networks
      if(mCp_props->contains("DirectedNetworks")){
	std::string netList = mCp_props->getProperty("DirectedNetworks");
	std::vector<std::string> directedNetworks;
	split(directedNetworks, netList, boost::is_any_of(" ,"), boost::token_compress_on);
	for(unsigned int i = 0; i < directedNetworks.size(); i++){
	  if(directedNetworks[i].compare("default") == 0) def.defineNetwork(true, &rlcm);
	  else                                	          def.defineNetwork(directedNetworks[i], true, &rlcm);
	}
      }

      if(mCp_props->contains("UndirectedNetworks")){
	std::string netList = mCp_props->getProperty("UndirectedNetworks");
	std::vector<std::string> undirectedNetworks;
	split(undirectedNetworks, netList, boost::is_any_of(" ,"), boost::token_compress_on);
	for(unsigned int i = 0; i < undirectedNetworks.size(); i++){
	  if(undirectedNetworks[i].compare("default") == 0) def.defineNetwork(false, &rlcm);
	  else                                              def.defineNetwork(undirectedNetworks[i], false, &rlcm);
	}
      }

      int processesPerX = strToInt(mCp_props->getProperty("proc.per.x"));
      int processesPerY = strToInt(mCp_props->getProperty("proc.per.y"));
      std::vector<int> procsPerDim;
      procsPerDim.push_back(processesPerX);
      procsPerDim.push_back(processesPerY);

      WorldCreator creator(comm);
      mCp_obs.reset(creator.createWorld<ObserverType, PatchType> (def, procsPerDim));
      mCp_obs->_setup(*mCp_props);
    }

    template<typename ObserverType, typename PatchType>
    void SimulationRunnerPlus<ObserverType,PatchType>::run() {
      RepastProcess::instance()->getScheduleRunner().run();
    }
  }

}

#endif /* SIMULATIONRUNNERPLUS_HPP_ */
