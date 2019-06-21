/**
 *  @file   LArReco/include/PandoraInterface.h
 * 
 *  @brief  Header file for PandoraInterface.
 * 
 *  $Log: $
 */
#ifndef PANDORA_INTERFACE_H
#define PANDORA_INTERFACE_H 1

#include "Pandora/PandoraInputTypes.h"

namespace pandora {class Pandora; class TiXmlElement;}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_reco
{

/**
 *  @brief  Parameters class
 */
class Parameters
{
public:
    /**
     *  @brief Default constructor
     */
    Parameters();

    std::string         m_settingsFile;                 ///< The path to the pandora settings file (mandatory parameter)
    std::string         m_eventFileNameList;            ///< Colon-separated list of file names to be processed

    int                 m_nEventsToProcess;             ///< The number of events to process (default all events in file)
    bool                m_shouldDisplayEventNumber;     ///< Whether event numbers should be displayed (default false)

    bool                m_shouldRunAllHitsCosmicReco;   ///< Whether to run all hits cosmic-ray reconstruction
    bool                m_shouldRunStitching;           ///< Whether to stitch cosmic-ray muons crossing between volumes
    bool                m_shouldRunCosmicHitRemoval;    ///< Whether to remove hits from tagged cosmic-rays
    bool                m_shouldRunSlicing;             ///< Whether to slice events into separate regions for processing
    bool                m_shouldRunNeutrinoRecoOption;  ///< Whether to run neutrino reconstruction for each slice
    bool                m_shouldRunCosmicRecoOption;    ///< Whether to run cosmic-ray reconstruction for each slice
    bool                m_shouldPerformSliceId;         ///< Whether to identify slices and select most appropriate pfos
    bool                m_printOverallRecoStatus;       ///< Whether to print current operation status messages

    float               m_wireAngleU;                   ///< Wire angle U
    float               m_wireAngleV;                   ///< Wire angle V
    float               m_wireAngleW;                   ///< Wire angle W
    float               m_wirePitchU;                   ///< TPC wire pitch U
    float               m_wirePitchV;                   ///< TPC wire pitch V
    float               m_wirePitchW;                   ///< TPC wire pitch W
    float               m_centerX;                      ///< TPC center x
    float               m_centerY;                      ///< TPC center y
    float               m_centerZ;                      ///< TPC center z
    float               m_widthX;                       ///< TPC width x
    float               m_widthY;                       ///< TPC width y
    float               m_widthZ;                       ///< TPC width z
    float               m_hitWidth;                     ///< Hit width in time

    pandora::InputInt   m_nEventsToSkip;                ///< The number of events to skip
};

/**
 *  @brief  ProtoHit class
 */
class ProtoHit
{
public:
    /**
     *  @brief Default constructor
     */
    ProtoHit();

    float               m_x;                            ///< Drift position
    float               m_z;                            ///< Wire number
    float               m_energy;                       ///< Energy
    pandora::HitType    m_hitType;                      ///< Hit type
    bool                m_deleteHit;                    ///< Has hit been merged and should original be deleted
};

typedef std::vector<ProtoHit> ProtoHitVector;
typedef std::map<pandora::HitType, float> HitTypeToFloatMap;

/**
 *  @brief  Create pandora instances
 *
 *  @param  parameters the parameters
 *  @param  pPrimaryPandora to receive the address of the primary pandora instance
 */
void CreatePandoraInstances(const Parameters &parameters, const pandora::Pandora *&pPrimaryPandora);

/**
 *  @brief  Process events using the supplied pandora instances
 *
 *  @param  parameters the application parameters
 *  @param  pPrimaryPandora the address of the primary pandora instance
 */
void ProcessEvents(const Parameters &parameters, const pandora::Pandora *const pPrimaryPandora);

/**
 *  @brief  Load geometry (hard coded for now)
 *
 *  @param  parameters the application parameters
 *  @param  pPrimaryPandora the address of the primary pandora instance
 */
void LoadGeometry(const Parameters &parameters, const pandora::Pandora *const pPrimaryPandora);

/**
 *  @brief  Load event information from xml file
 *
 *  @param  parameters the application parameters
 *  @param  pPrimaryPandora the address of the primary pandora instance
 *  @param  nEvents event number to load
 *  @param  pTiXmlElement xml element
 */
void LoadEvent(const Parameters &parameters, const pandora::Pandora *const pPrimaryPandora, const int nEvents, TiXmlElement *pTiXmlElement);

/**
 *  @brief  Load hit from xml file
 *
 *  @param  pTiXmlElement xml element
 *  @param  protoHitVectorU U view proto hits
 *  @param  protoHitVectorV V view proto hits
 *  @param  protoHitVectorW W view proto hits
 */
void LoadCell(TiXmlElement *pTiXmlElement, ProtoHitVector &protoHitVectorU, ProtoHitVector &protoHitVectorV, ProtoHitVector &protoHitVectorW);

/**
 *  @brief  Downsample hits
 *
 *  @param  parameters the application parameters
 *  @param  protoHitVector vector of protoHits
 */
void DownsampleHits(const Parameters &inputParameters, ProtoHitVector &protoHitVector);

/**
 *  @brief  Search and return hits to merge
 *
 *  @param  parameters the application parameters
 *  @param  protoHitVector vector of protoHits
 *  @param  protoHit1 merge candidate one
 *  @param  protoHit2 merge candidate two
 *
 *  @return is a merge is present
 */
bool IdentifyMerge(const Parameters &inputParameters, ProtoHitVector &protoHitVector, ProtoHit &protoHit1, ProtoHit &protoHit2);

/**
 *  @brief  Convert the YZ position to a U
 *
 *  @param  y position
 *  @param  z position
 *  @param  parameters the application parameters
 *
 *  @return the u position
 */
float YZtoU(const float y, const float z, const Parameters &parameters);

/**
 *  @brief  Convert the YZ position to a V
 *
 *  @param  y position
 *  @param  z position
 *  @param  parameters the application parameters
 *
 *  @return the v position
 */
float YZtoV(const float y, const float z, const Parameters &parameters);

/**
 *  @brief  Sort ProtoHits by position
 *
 *  @param  protoHit1 first hit
 *  @param  protoHit2 second hit
 *
 *  @return is protoHit1 sorted above protoHit2
 */
bool SortProtoHits(const ProtoHit &protoHit1, const ProtoHit &protoHit2);

/**
 *  @brief  Parse the command line arguments, setting the application parameters
 *
 *  @param  argc argument count
 *  @param  argv argument vector
 *  @param  parameters to receive the application parameters
 *
 *  @return success
 */
bool ParseCommandLine(int argc, char *argv[], Parameters &parameters);

/**
 *  @brief  Load xml element
 *
 *  @param  pHeadTiXmlElement pointer to xml element
 *  @param  value to set
 *  @param  name of element
 */
void LoadXmlElement(const pandora::TiXmlElement *pHeadTiXmlElement, float &value, const std::string &name);

/**
 *  @brief  Print the list of configurable options
 *
 *  @return false, to force abort
 */
bool PrintOptions();

/**
 *  @brief  Process the provided reco option string to perform high-level steering
 *
 *  @param  recoOption the reco option string
 *  @param  parameters to receive the application parameters
 *
 *  @return success
 */
bool ProcessRecoOption(const std::string &recoOption, Parameters &parameters);

/**
 *  @brief  Process list of external, commandline parameters to be passed to specific algorithms
 *
 *  @param  parameters the parameters
 *  @param  pPandora the address of the pandora instance
 */
void ProcessExternalParameters(const Parameters &parameters, const pandora::Pandora *const pPandora);

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline Parameters::Parameters() :
    m_settingsFile(""),
    m_eventFileNameList(""),
    m_nEventsToProcess(-1),
    m_shouldDisplayEventNumber(false),
    m_shouldRunAllHitsCosmicReco(true),
    m_shouldRunStitching(true),
    m_shouldRunCosmicHitRemoval(true),
    m_shouldRunSlicing(true),
    m_shouldRunNeutrinoRecoOption(true),
    m_shouldRunCosmicRecoOption(true),
    m_shouldPerformSliceId(true),
    m_printOverallRecoStatus(false),
    m_wireAngleU(std::numeric_limits<float>::max()),
    m_wireAngleV(std::numeric_limits<float>::max()),
    m_wireAngleW(std::numeric_limits<float>::max()),
    m_wirePitchU(std::numeric_limits<float>::max()),
    m_wirePitchV(std::numeric_limits<float>::max()),
    m_wirePitchW(std::numeric_limits<float>::max()),
    m_centerX(std::numeric_limits<float>::max()),
    m_centerY(std::numeric_limits<float>::max()),
    m_centerZ(std::numeric_limits<float>::max()),
    m_widthX(std::numeric_limits<float>::max()),
    m_widthY(std::numeric_limits<float>::max()),
    m_widthZ(std::numeric_limits<float>::max()),
    m_hitWidth(std::numeric_limits<float>::max())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline ProtoHit::ProtoHit() :
    m_x(std::numeric_limits<float>::max()),
    m_z(std::numeric_limits<float>::max()),
    m_energy(std::numeric_limits<float>::max()),
    m_hitType(pandora::TPC_3D),
    m_deleteHit(false)
{
}

} // namespace lar_reco

#endif // #ifndef PANDORA_INTERFACE_H
