/*
 * Repast for High Performance Computing (Repast HPC)
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
 *  Created on: Jun 14, 2011
 *      Author: John T. Murphy
 */

#ifndef AGENTIMPORTEREXPORTER_H_
#define AGENTIMPORTEREXPORTER_H_

#include <map>
#include <set>
#include <list>

#include "AgentRequest.h"
#include "AgentId.h"
#include "AgentStatus.h"
#include "mpi_constants.h"

#define SHARE_AGENTS_BY_SET

// If you wish to allow a single cancellation to be used
// to cancel all of an agent's requests, define this
//#define ALLOW_FULL_AGENT_REQUEST_CANCELLATION


/* Preprocessor configuration */

#ifndef DEFAULT_IMPORTER_EXPORTER
  #define DEFAULT_IMPORTER_EXPORTER 1
#endif

#if DEFAULT_IMPORTER_EXPORTER < 1 || DEFAULT_IMPORTER_EXPORTER > 5
  #error "Invalid value used for default Importer_Exporter"
#endif

// If you are not permitting agent request by set,
// then only one of the importer_exporters should be compiled
// (The one chosen as the default)
#ifndef SHARE_AGENTS_BY_SET

  #if DEFAULT_IMPORTER_EXPORTER != 1 && !defined OMIT_IMPORTER_EXPORTER_COUNT_LIST
  #define OMIT_IMPORTER_EXPORTER_COUNT_LIST
  #endif

  #if DEFAULT_IMPORTER_EXPORTER != 2 && !defined OMIT_IMPORTER_EXPORTER_COUNT_SET
  #define OMIT_IMPORTER_EXPORTER_COUNT_SET
  #endif

  #if DEFAULT_IMPORTER_EXPORTER != 3 && !defined OMIT_IMPORTER_EXPORTER_LIST
  #define OMIT_IMPORTER_EXPORTER_LIST
  #endif

  #if DEFAULT_IMPORTER_EXPORTER != 4 && !defined OMIT_IMPORTER_EXPORTER_SET
  #define OMIT_IMPORTER_EXPORTER_SET
  #endif

  #if DEFAULT_IMPORTER_EXPORTER != 5 && !defined OMIT_IMPORTER_EXPORTER_MAP_int
  #define OMIT_IMPORTER_EXPORTER_MAP_int
  #endif

#endif



// Cannot set both 'default' and 'omit'; setting one
// as default overrides the omit
#if DEFAULT_IMPORTER_EXPORTER == 1
  #define DEFAULT_IMPORTER_EXPORTER_CLASS ImporterExporter_COUNT_LIST
  #define DEFAULT_ENUM_SYMBOL COUNT_LIST
  #undef OMIT_IMPORTER_EXPORTER_COUNT_LIST
#endif

#if DEFAULT_IMPORTER_EXPORTER == 2
  #define DEFAULT_IMPORTER_EXPORTER_CLASS ImporterExporter_COUNT_SET
  #define DEFAULT_ENUM_SYMBOL COUNT_SET
  #undef OMIT_IMPORTER_EXPORTER_COUNT_SET
#endif

#if DEFAULT_IMPORTER_EXPORTER == 3
  #define DEFAULT_IMPORTER_EXPORTER_CLASS ImporterExporter_LIST
  #define DEFAULT_ENUM_SYMBOL LIST
  #undef OMIT_IMPORTER_EXPORTER_LIST
#endif

#if DEFAULT_IMPORTER_EXPORTER == 4
  #define DEFAULT_IMPORTER_EXPORTER_CLASS ImporterExporter_SET
  #define DEFAULT_ENUM_SYMBOL SET
  #undef OMIT_IMPORTER_EXPORTER_SET
#endif

#if DEFAULT_IMPORTER_EXPORTER == 5
  #define DEFAULT_IMPORTER_EXPORTER_CLASS ImporterExporter_Map_int
  #define DEFAULT_ENUM_SYMBOL MAP_int
  #undef OMIT_IMPORTER_EXPORTER_MAP_int
#endif



#ifdef SHARE_AGENTS_BY_SET
// This enum is useful when requests by set are allowed
// Result of this construction is enum AGENT_IMPORTER_EXPORTER_TYPE{ Default[, non-default1][, non-default2] }, with non-compiled classes omitted
enum AGENT_IMPORTER_EXPORTER_TYPE{ DEFAULT_ENUM_SYMBOL
#if !defined OMIT_IMPORTER_EXPORTER_COUNT_LIST && DEFAULT_IMPORTER_EXPORTER != 1
, COUNT_LIST
#endif
#if !defined OMIT_IMPORTER_EXPORTER_COUNT_SET && DEFAULT_IMPORTER_EXPORTER != 2
, COUNT_SET
#endif
#if !defined OMIT_IMPORTER_EXPORTER_LIST && DEFAULT_IMPORTER_EXPORTER != 3
, LIST
#endif
#if !defined OMIT_IMPORTER_EXPORTER_SET && DEFAULT_IMPORTER_EXPORTER != 4
, SET
#endif
#if !defined OMIT_IMPORTER_EXPORTER_MAP_int && DEFAULT_IMPORTER_EXPORTER != 5
, MAP_int
#endif
};

#define DEFAULT_AGENT_REQUEST_SET "__Default_Agent_Request_Set__"

#endif

/* END PREPROCESSOR */


#define REQUEST_AGENTS_ALL "AGENTS_FROM_ALL_IMPORTER_EXPORTERS"


namespace repast{


#ifndef SHARE_AGENTS_BY_SET
typedef std::map<int, AgentRequest*>  AgentExporterInfo;
#else
typedef std::map<std::string, std::map<int, AgentRequest*>* > AgentExporterInfo;
#endif

/**
 * Data structure for exporter data that is to be sent
 * to other processes when the agents being exported are
 * moved. Note that the internal data structure is protected;
 * classes can use this data without knowing its actual internal
 * structure.
 */
class AgentExporterData{
private:
  AgentExporterInfo empty; // Used as dummy for sending empty data
  std::map<int, AgentExporterInfo* > data;

#ifdef SHARE_AGENTS_BY_SET
  std::string selectedSet;
#endif

public:
  AgentExporterData();
  ~AgentExporterData();

  /**
   * Adds an agent ID to this list of data that is being exported
   * to a specific processor (destProc), so that the agent's
   * information will be exported to another processor (sourceProc).
   */
  void addData(const AgentId& id, const int destProc, const int sourceProc,
      const int numberOfCopies = 1);

  /**
   * Gets the packaged set of information to be sent to a specific
   * processor
   */
  AgentExporterInfo* dataForProc(int destProc);

  /**
   * Clears this data structure
   */
  void clear();

  /**
   * Remove all the data for a specific agent; useful
   * when the agent is removed
   */
  void removeAllDataForAgent(AgentId& id);


#ifdef SHARE_AGENTS_BY_SET
  /**
   * Specifies that add and retrieve actions are to be
   * performed on the subset of data identified by the given set
   * name. (Does not affect 'clear' or 'removeAllDataForAgent')
   */
  void selectSet(std::string setName);
#endif

};




/* Importers */

#ifdef ALLOW_FULL_AGENT_REQUEST_CANCELLATION
/**
 * Marks an ID as a 'full cancellation,' meaning
 * all previous requests for this agent will be
 * cancelled.
 */
void markIdAsFullCancellation(AgentId& id){
  id.currentRank(((-1) * id.currentRank()) - 1);
}

/**
  * Checks an AgentId to see if it is marked as a 'full cancellation';
  * if it is, MODIFIES IT back to its original, unmarked
  * state and returns true
  *
  * @param id the AgentId to be checked; modified back to
  * a normal AgentId if found to be a full cancellation
  * @return true if the AgentId was in a state that marked it
  * as a full cancellation, false otherwise
  */
bool checkForFullCancellation(AgentId& id){
  if(id.currentRank() < 0){
    id.currentRank((id.currentRank() + 1) * (-1));
    return true;
  }
  else return false;
}
#endif



/* IMPORTER */

/**
 * This class manages importing agent information; primarily
 * this means constructing the appropriate mpi receives when
 * agent information is to be exchanged. However, this class
 * can also define specific semantics that can apply to
 * agent requests- what to do in the case that an agent
 * is requested twice, for example.
 */
class AbstractImporter{
protected:
  // Each child implementation is required to maintain this list
  std::set<int> exportingProcesses;

public:

  AbstractImporter();
  virtual ~AbstractImporter();

  /**
   * Gets a const reference to the set of ints representing the processes that are
   * sending this process agent information
   */
  virtual const std::set<int>& getExportingProcesses(){ return exportingProcesses; }

  /**
   * Given an agent request (including requests for agents on multiple other
   * processes), makes a record of the agents that are being requested by this
   * process and will therefore be received from other processes. The record must
   * at a minimum indicate which other processes are sending agent information,
   * but may include other information, such as how many times a particular
   * agent has been requested.
   */
  virtual void registerOutgoingRequests(AgentRequest &req) = 0;

  /**
   * Notifies this importer that the agent that it (presumably) has been importing
   * has been removed from the simulation on its home process, and the information
   * for that agent will no longer be sent
   */
  virtual void importedAgentIsRemoved(const AgentId& id) = 0;

  /**
   * Notifies this importer that the agent that it (presumably) has been importing
   * from another process has been moved; its information will now be received
   * from its new home process (unless the agent was moved to this process)
   */
  virtual void importedAgentIsMoved(const AgentId& id, int newProcess) = 0;

  /**
   * Some semantic sugar; operationally this is the same as 'importedAgentIsMoved'
   */
  inline void importedAgentIsNowLocal(const AgentId& id){ importedAgentIsRemoved(id); }

 };


#if !defined OMIT_IMPORTER_EXPORTER_COUNT_LIST || \
    !defined OMIT_IMPORTER_EXPORTER_COUNT_SET
class Importer_COUNT: public AbstractImporter{

private:
 std::map<int, int> sources;
 void incrementCount(int sourceProcess);
 void decrementCount(int sourceProcess);

public:
  Importer_COUNT();
  virtual ~Importer_COUNT();

  virtual void registerOutgoingRequests(AgentRequest &req);
  virtual void importedAgentIsRemoved(const AgentId& id);
  virtual void importedAgentIsMoved(const AgentId& id, int newProcess);

};
#endif

#ifndef OMIT_IMPORTER_EXPORTER_LIST
class Importer_LIST: public AbstractImporter{

private:
  std::map<int, std::list<AgentId>* > sources;

  // Inline helper functions to manage the internal recordkeeping

  // Gets the record associated with the specified other process;
  // if none exists, one is created. 'exportingProcesses' is updated
  inline std::list<AgentId>* getRecord(int rank);

  // Removes all the entries for a given agent from the specified record;
  // Returns the number of entries removed
  // Also: Invokes 'checkRecord'
  inline int removeAll(const AgentId& id, std::list<AgentId>* record, int rank);

  // Removes one entry for a given ID (passes through to removeId by position)
  inline void removeID(AgentId& id, std::list<AgentId>* record, int rank);

  // Removes the entry at the specified position in the record list.
  // Also: Invokes 'checkRecord'
  inline void removeID(std::list<AgentId>::iterator idPosition, std::list<AgentId>* record, int rank);

  // Checks the length of this record; if the record is now empty,
  // it is removed from the map (and exportingProcesses is updated)
  inline void checkRecord(std::list<AgentId>* record, int rank);

public:
  Importer_LIST();
  ~Importer_LIST();

  virtual void registerOutgoingRequests(AgentRequest& req);
  virtual void importedAgentIsRemoved(const AgentId& id);
  virtual void importedAgentIsMoved(const AgentId& id, int newProcess);

};
#endif


#ifndef OMIT_IMPORTER_EXPORTER_SET
class Importer_SET: public AbstractImporter{
private:
  std::map<int, std::set<AgentId>* > sources;

  // Helper functions to manage internal bookkeeping

  // Gets the record associated with the specified other
  // rank. If none exists, one is created (and exportingProcesses
  // is updated).
  inline std::set<AgentId>* getRecord(int rank);

  // Removes an ID from the specified set.
  // If this results in an actual removal, will return 1; if not,
  // will return zero. AgentID will not be in map after operation.
  // Also: Invokes checkRecord
  inline int removeID(const AgentId& id, std::set<AgentId>* record, int rank);

  // Checks the size of the record passed; if zero, deletes it from
  // the map and updates exportingProcesses
  inline void checkRecord(std::set<AgentId>* record, int rank);

public:
  Importer_SET();
  ~Importer_SET();

  virtual void registerOutgoingRequests(AgentRequest& req);
  virtual void importedAgentIsRemoved(const AgentId& id);
  virtual void importedAgentIsMoved(const AgentId& id, int newProcess);

};
#endif


#ifndef OMIT_IMPORTER_EXPORTER_MAP_int
class Importer_MAP_int: public AbstractImporter{
private:
  std::map<int, std::map<AgentId, int>* > sources;

  // Helper functions to manage internal bookkeeping

  // Gets the record associated with the specified rank
  // If one does not exist, it is created, and exportingProcesses
  // is updated
  std::map<AgentId, int>* getRecord(int rank);

  // Removes the entry for the specified agent.
  // Returns the number associated with that agent (essentially the number of 'entries' removed)
  // Also: invokes checkRecord
  inline int removeAll(const AgentId& id, std::map<AgentId, int>* record, int rank);

  // Removes a single entry for the specified ID (decrements its count)
  // If the entry does not exist returns 0, otherwise
  // returns one (representing one entry removed)
  // Also: invokes checkRecord
  inline int removeID(const AgentId& id, std::map<AgentId, int>* record, int rank);

  // Removes the entire entry for a given agent, by position. Does not return
  // any indication of success of operation
  // Invokes checkRecord
  inline void removeID(std::map<AgentId, int>::iterator idPosition, std::map<AgentId, int>* record, int rank);

  // Checks if the record size is zero; if so, deletes it and notifies exportingProcesses
  inline void checkRecord(std::map<AgentId, int>* record, int rank);


public:
  Importer_MAP_int();
  ~Importer_MAP_int();

  virtual void registerOutgoingRequests(AgentRequest& req);
  virtual void importedAgentIsRemoved(const AgentId& id);
  virtual void importedAgentIsMoved(const AgentId& id, int newProcess);

};
#endif










/* Exporters */

static std::vector<AgentStatus>    emptyStatus;

/**
 * The 'Exporter' class is responsible for keeping a list
 * of the agents that have been requested by other processes
 * and whose data is to be sent to them when agents' states
 * are synchronized, and for packaging and sending
 * that data during synchronization. It is also responsible
 * for exchanging this 'export' information when any of the
 * agents that it is exporting are being moved to other
 * processes; when an agent moves, its new home process must
 * be able to assume the same export duties that its original
 * process was performing.
 */
class AbstractExporter{
public:
  typedef std::map<int, std::set<AgentStatus> > StatusMap;

private:
  StatusMap*                  outgoingStatusChangesDeletePtr;
  AgentExporterData*          outgoingAgentExporterInformationDeletePtr;


protected:
  StatusMap*                  outgoingStatusChanges;
  AgentExporterData*          outgoingAgentExporterInformation;

  std::set<int>               processesExportedTo;
  std::map<int, AgentRequest> exportedMap;

public:

  AbstractExporter();
  virtual ~AbstractExporter();

  /**
   * Makes a record of the data receives (in the form of a vector of AgentRequests)
   * so that the agents' data can be sent to the requesting processes.
   */
  virtual void registerIncomingRequests(std::vector<AgentRequest>& requests) = 0;

  /**
   * The set of information received here comprises the information that some
   * other process was using to export information about agents that are
   * now being moved to this process. This method takes that information
   * and incorporates it into this exporter, so that this exporter can now
   * export the agents' information to the processes that have requested it.
   */
  virtual void incorporateAgentExporterInfo(std::map<int, AgentRequest*> info);


  /**
   * 1) Removes the agent export information from this process
   * 2) Updates the outgoing status change buffer to include
   *      the status change for this agent to all procs
   *      to which this agent was being exported (except
   *      if one of these was the proc to which the agent
   *      is now moving; this is omitted)
   */
  virtual void agentRemoved(const AgentId& id);

  /**
   * 1) Removes the agent export information from this process
   * 2) Places a copy of the agent export information into the
   *      outgoing buffer
   * 3) Updates the outgoing status change buffer to include
   *      the status change for this agent to all procs
   *      to which this agent was being exported (except
   *      if one of these was the proc to which the agent
   *      is now moving; this is omitted)
   */
  virtual void agentMoved(const AgentId& id, int process);

  /**
   * Gets the list of processes this exporter is sending information
   * to.
   */
  virtual const std::set<int>& getProcessesExportedTo();

  /**
   * Gets the export information that has been placed into the 'outgoing agent
   * export information' buffer because agents that were being exported are being
   * sent to a new process, for the specified process.
   */
  AgentExporterInfo* getAgentExportInfo(int destProc);

  /**
   * Gets the set of status changes for the exported agents.
   */
  const StatusMap* getOutgoingStatusChanges();

  /**
   * Clears the outgoing agent export information buffer; should be
   * called after the information is sent.
   */
  void clearAgentExportInfo();

  /**
   * Clears the outgoing status information buffer; should
   * be called after the information is sent
   */
  void clearStatusMap();

  /**
   * Gets the list of agents being exported by this exported,
   * as a map by ints representing the processes to which
   * information will be sent.
   */
  virtual const std::map<int, AgentRequest>& getAgentsToExport();

#ifdef SHARE_AGENTS_BY_SET
public:
  AbstractExporter(StatusMap* outgoingStatusMap, AgentExporterData* outgoingAgentExporterInfo);
#endif
};


#if !defined OMIT_IMPORTER_EXPORTER_COUNT_LIST || \
    !defined OMIT_IMPORTER_EXPORTER_LIST       || \
    !defined OMIT_IMPORTER_EXPORTER_SET        || \
    !defined OMIT_IMPORTER_EXPORTER_MAP_int
class Exporter_LIST: public AbstractExporter{
private:

public:
  Exporter_LIST();
  virtual ~Exporter_LIST();

  virtual void registerIncomingRequests(std::vector<AgentRequest>& requests);

#ifdef SHARE_AGENTS_BY_SET
public:
  Exporter_LIST(StatusMap* outgoingStatusMap, AgentExporterData* outgoingAgentExporterInfo);
#endif
};
#endif


#ifndef OMIT_IMPORTER_EXPORTER_COUNT_SET
class Exporter_SET: public AbstractExporter{

public:
  Exporter_SET();
  virtual ~Exporter_SET();

  virtual void registerIncomingRequests(std::vector<AgentRequest>& requests);

#ifdef SHARE_AGENTS_BY_SET
public:
  Exporter_SET(StatusMap* outgoingStatusMap, AgentExporterData* outgoingAgentExporterInfo);
#endif
};
#endif


/* Importer and Exporter */
class AbstractImporterExporter{

protected:
  AbstractImporter* importer;
  AbstractExporter* exporter;

public:
  AbstractImporterExporter(AbstractImporter* i, AbstractExporter* e);
  virtual ~AbstractImporterExporter();

  // Pass-through methods for all importer and exporter public methods
  virtual const std::set<int>& getExportingProcesses(){                                            return importer->getExportingProcesses();               }
  virtual       void           registerOutgoingRequests(AgentRequest &req){                               importer->registerOutgoingRequests(req);         }
  virtual       void           importedAgentIsRemoved(const AgentId& id){                                 importer->importedAgentIsRemoved(id);            }
  virtual       void           importedAgentIsMoved(const AgentId& id, int newProcess){                   importer->importedAgentIsMoved(id, newProcess);  }

  inline void importedAgentIsNowLocal(const AgentId& id){                                                 importer->importedAgentIsNowLocal(id);           }

  virtual const AbstractExporter::StatusMap* getOutgoingStatusChanges();

  virtual const std::set<int>& getProcessesExportedTo(){                                           return exporter->getProcessesExportedTo();              }
  virtual       void           registerIncomingRequests(std::vector<AgentRequest>& requests){             exporter->registerIncomingRequests(requests);    }
  virtual       void           agentRemoved(const AgentId& id){                                           exporter->agentRemoved(id);                      }
  virtual       void           agentMoved(const AgentId& id, int process){                                exporter->agentMoved(id, process);               }
  virtual       void           incorporateAgentExporterInfo(std::map<int, AgentRequest* > info){          exporter->incorporateAgentExporterInfo(info);    }
  virtual       void           clearStatusMap(){                                                          exporter->clearStatusMap();                      }
  virtual AgentExporterInfo*   getAgentExportInfo(int destProc){                                   return exporter->getAgentExportInfo(destProc);          }
  virtual       void           clearAgentExportInfo(){                                                    exporter->clearAgentExportInfo();                }
  virtual const std::map<int, AgentRequest>& getAgentsToExport(){                                  return exporter->getAgentsToExport();                   }

  /**
   * Exchanges the contents of the 'statusMap' with the destination processes, updating
   * the status (moved or removed) for all agents being exported. Returns this information
   * in the statusUpdates vector.
   */
  virtual void exchangeAgentStatusUpdates(boost::mpi::communicator world, std::vector<std::vector<AgentStatus>* >& statusUpdates);

};

/* Normal variants, with semantics defined by which importer/exporter combination is used */

#ifndef OMIT_IMPORTER_EXPORTER_COUNT_LIST
class ImporterExporter_COUNT_LIST: public AbstractImporterExporter{
public:
  ImporterExporter_COUNT_LIST();

#ifdef SHARE_AGENTS_BY_SET
  ImporterExporter_COUNT_LIST(AbstractExporter::StatusMap* outgoingStatusMap, AgentExporterData* outgoingAgentExporterInfo);
#endif

  virtual ~ImporterExporter_COUNT_LIST();

};
#endif

#ifndef OMIT_IMPORTER_EXPORTER_COUNT_SET
class ImporterExporter_COUNT_SET: public AbstractImporterExporter{
public:
  ImporterExporter_COUNT_SET();

#ifdef SHARE_AGENTS_BY_SET
  ImporterExporter_COUNT_SET(AbstractExporter::StatusMap* outgoingStatusMap, AgentExporterData* outgoingAgentExporterInfo);
#endif

  virtual ~ImporterExporter_COUNT_SET();
};
#endif

#ifndef OMIT_IMPORTER_EXPORTER_LIST
class ImporterExporter_LIST: public AbstractImporterExporter{
public:
  ImporterExporter_LIST();

#ifdef SHARE_AGENTS_BY_SET
  ImporterExporter_LIST(AbstractExporter::StatusMap* outgoingStatusMap, AgentExporterData* outgoingAgentExporterInfo);
#endif

  virtual ~ImporterExporter_LIST();
};
#endif

#ifndef OMIT_IMPORTER_EXPORTER_SET
class ImporterExporter_SET: public AbstractImporterExporter{
public:
  ImporterExporter_SET();

#ifdef SHARE_AGENTS_BY_SET
  ImporterExporter_SET(AbstractExporter::StatusMap* outgoingStatusMap, AgentExporterData* outgoingAgentExporterInfo);
#endif

  virtual ~ImporterExporter_SET();
};
#endif

#ifndef OMIT_IMPORTER_EXPORTER_MAP_int
class ImporterExporter_MAP_int: public AbstractImporterExporter{
public:
  ImporterExporter_MAP_int();

#ifdef SHARE_AGENTS_BY_SET
  ImporterExporter_MAP_int(AbstractExporter::StatusMap* outgoingStatusMap, AgentExporterData* outgoingAgentExporterInfo);
#endif

  virtual ~ImporterExporter_MAP_int();
};
#endif


/* "BY SET" variant; allows multiple sets of shared agents to be managed independently */

#ifdef SHARE_AGENTS_BY_SET
class ImporterExporter_BY_SET: public AbstractImporterExporter{

private:

  AbstractExporter::StatusMap*                      outgoingStatusChanges;
  AgentExporterData*                                outgoingAgentExporterInformation;
  std::map<std::string, AbstractImporterExporter*>  importersExportersMap;
  std::set<int>                                     exportingProcesses;
  std::set<int>                                     processesExportedTo;
  std::map<int, AgentRequest>                       exportedMap;

  AbstractImporterExporter* getSet(std::string setname,
      AGENT_IMPORTER_EXPORTER_TYPE setType = DEFAULT_ENUM_SYMBOL);


  inline void rebuildExportingProcesses();

  inline void rebuildProcessesExportedTo();

  inline void rebuildExportedMap();

public:
  ImporterExporter_BY_SET();

  virtual ~ImporterExporter_BY_SET();


  /* Importer-related functions */

  virtual const std::set<int>& getExportingProcesses();
          const std::set<int>& getExportingProcesses(std::string setName);


  virtual void registerOutgoingRequests(AgentRequest& request);
          void registerOutgoingRequests(AgentRequest& request, std::string setName,
                     AGENT_IMPORTER_EXPORTER_TYPE setType = DEFAULT_ENUM_SYMBOL);

  virtual void importedAgentIsRemoved(const AgentId& id);

  virtual void importedAgentIsMoved(const AgentId& id, int newProcess);


  /* Exporter-related functions */

  virtual const AbstractExporter::StatusMap* getOutgoingStatusChanges();

  virtual const std::set<int>& getProcessesExportedTo();
          const std::set<int>& getProcessesExportedTo(std::string setName);

  virtual void registerIncomingRequests(std::vector<AgentRequest>& requests);
          void registerIncomingRequests(std::vector<AgentRequest>& requests, std::string setName);

  virtual void agentRemoved(const AgentId& id);
  virtual void agentMoved(const AgentId& id, int newProcess);

  // This should NEVER be called; it is inappropriate to use this when 'by set' is being used
  virtual void incorporateAgentExporterInfo(std::map<int, AgentRequest* > info);
          void incorporateAgentExporterInfo(std::map<std::string, std::map<int, AgentRequest*>*>  info);

  virtual void clearStatusMap();

  virtual AgentExporterInfo* getAgentExportInfo(int destProc);

  virtual void clearAgentExportInfo();

  virtual const std::map<int, AgentRequest>& getAgentsToExport();
          const std::map<int, AgentRequest>& getAgentsToExport(std::string setName);


  void dropSet(std::string setName){
    importersExportersMap.erase(setName); // Fails silently if the specified set is not present
  }


};
#endif

}
#endif /* AGENTIMPORTEREXPORTER_H_ */
