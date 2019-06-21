/**
 *  @file   LArRecoMP/test/PandoraInterface.cc
 *
 *  @brief  Implementation of the lar reco mp application
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"
#include "Helpers/XmlHelper.h"
#include "Xml/tinyxml.h"

#include "larpandoracontent/LArContent.h"
#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"
#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArPersistency/EventReadingAlgorithm.h"
#include "larpandoracontent/LArPlugins/LArPseudoLayerPlugin.h"
#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"

#ifdef LIBTORCH_DL
#include "larpandoradlcontent/LArDLContent.h"
#endif

#include "PandoraInterface.h"

#ifdef MONITORING
#include "TApplication.h"

#include "TChain.h"
#endif

#include <getopt.h>
#include <iostream>
#include <string>

using namespace pandora;
using namespace lar_reco;

int main(int argc, char *argv[])
{
    int errorNo(0);
    const Pandora *pPrimaryPandora(nullptr);

    try
    {
        Parameters parameters;

        if (!ParseCommandLine(argc, argv, parameters))
            return 1;

#ifdef MONITORING
        TApplication *pTApplication = new TApplication("LArReco", &argc, argv);
        pTApplication->SetReturnFromRun(kTRUE);
#endif
        CreatePandoraInstances(parameters, pPrimaryPandora);

        if (!pPrimaryPandora)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        ProcessEvents(parameters, pPrimaryPandora);
    }
    catch (const StatusCodeException &statusCodeException)
    {
        std::cerr << "Pandora StatusCodeException: " << statusCodeException.ToString() << statusCodeException.GetBackTrace() << std::endl;
        errorNo = 1;
    }
    catch (const StopProcessingException &stopProcessingException)
    {
        // Exit gracefully
        std::cout << stopProcessingException.GetDescription() << std::endl;
        errorNo = 0;
    }
    catch (...)
    {
        std::cerr << "Unknown exception: " << std::endl;
        errorNo = 1;
    }

    MultiPandoraApi::DeletePandoraInstances(pPrimaryPandora);
    return errorNo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_reco
{

void CreatePandoraInstances(const Parameters &parameters, const Pandora *&pPrimaryPandora)
{
    pPrimaryPandora = new Pandora();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterAlgorithms(*pPrimaryPandora));
#ifdef LIBTORCH_DL
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLContent::RegisterAlgorithms(*pPrimaryPandora));
#endif
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterBasicPlugins(*pPrimaryPandora));

    if (!pPrimaryPandora)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    MultiPandoraApi::AddPrimaryPandoraInstance(pPrimaryPandora);

    ProcessExternalParameters(parameters, pPrimaryPandora);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(*pPrimaryPandora, new lar_content::LArPseudoLayerPlugin));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetLArTransformationPlugin(*pPrimaryPandora, new lar_content::LArRotationalTransformationPlugin));
    LoadGeometry(parameters, pPrimaryPandora);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*pPrimaryPandora, parameters.m_settingsFile));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProcessEvents(const Parameters &parameters, const Pandora *const pPrimaryPandora)
{
    int nEvents(0);

    // Load Input File
    TiXmlDocument *pTiXmlDocument = new TiXmlDocument();

    if (!pTiXmlDocument->LoadFile(parameters.m_eventFileNameList.c_str()))
        std::cerr << pTiXmlDocument->ErrorDesc() << std::endl;

    TiXmlElement *pRunTiXmlElement = pTiXmlDocument->FirstChildElement();
    TiXmlElement *pEventTiXmlElement = pRunTiXmlElement->FirstChildElement();

    while ((nEvents++ < parameters.m_nEventsToProcess) || (0 > parameters.m_nEventsToProcess))
    {
        if (!pEventTiXmlElement)
        {
            pTiXmlDocument->Clear();
            delete pTiXmlDocument;
            delete pRunTiXmlElement;
            delete pEventTiXmlElement;
            break;
        }

        if (parameters.m_shouldDisplayEventNumber)
            std::cout << std::endl << "   PROCESSING EVENT: " << (nEvents - 1) << std::endl << std::endl;

        LoadEvent(parameters, pPrimaryPandora, pEventTiXmlElement);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*pPrimaryPandora));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*pPrimaryPandora));

        pEventTiXmlElement->NextSiblingElement();
        pEventTiXmlElement->NextSiblingElement();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LoadGeometry(const Parameters &inputParameters, const Pandora *const pPrimaryPandora)
{
    try
    {
        PandoraApi::Geometry::LArTPC::Parameters parameters;
        parameters.m_larTPCVolumeId = 0;
        parameters.m_centerX = inputParameters.m_centerX;
        parameters.m_centerY = inputParameters.m_centerY;
        parameters.m_centerZ = inputParameters.m_centerZ;
        parameters.m_widthX = inputParameters.m_widthX;
        parameters.m_widthY = inputParameters.m_widthY;
        parameters.m_widthZ = inputParameters.m_widthZ;
        parameters.m_wirePitchU = inputParameters.m_wirePitchU;
        parameters.m_wirePitchV = inputParameters.m_wirePitchV;
        parameters.m_wirePitchW = inputParameters.m_wirePitchW;
        parameters.m_wireAngleU = inputParameters.m_wireAngleU;
        parameters.m_wireAngleV = inputParameters.m_wireAngleV;
        parameters.m_wireAngleW = inputParameters.m_wireAngleW;
        parameters.m_sigmaUVW = 1.51300001144;
        parameters.m_isDriftInPositiveX = true;
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LArTPC::Create(*pPrimaryPandora, parameters));
    }
    catch (...)
    {
        std::cout << "Cannot load geometry" << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LoadEvent(const Parameters &inputParameters, const pandora::Pandora *const pPrimaryPandora, TiXmlElement *pTiXmlElement)
{
    ProtoHitVector protoHitVectorU, protoHitVectorV, protoHitVectorW;

    for (TiXmlElement *pSubTiXmlElement = pTiXmlElement->FirstChildElement(); pSubTiXmlElement != NULL; pSubTiXmlElement = pSubTiXmlElement->NextSiblingElement())
    {
        const std::string componentName(pSubTiXmlElement->ValueStr());

        if (componentName == "Cell")
        {
            LoadCell(inputParameters, pSubTiXmlElement, protoHitVectorU, protoHitVectorV, protoHitVectorW);
        }
        else if (componentName == "MCParticle")
        {
            LoadMCParticle(pSubTiXmlElement, pPrimaryPandora);
        }
    }

    DownsampleHits(inputParameters, protoHitVectorU);
    DownsampleHits(inputParameters, protoHitVectorV);
    DownsampleHits(inputParameters, protoHitVectorW);

    ProtoHitVector protoHitVector;
    protoHitVector.insert(protoHitVector.end(), protoHitVectorU.begin(), protoHitVectorU.end());
    protoHitVector.insert(protoHitVector.end(), protoHitVectorV.begin(), protoHitVectorV.end());
    protoHitVector.insert(protoHitVector.end(), protoHitVectorW.begin(), protoHitVectorW.end());

    HitTypeToFloatMap hitTypeToThickness;
    hitTypeToThickness.insert(HitTypeToFloatMap::value_type(pandora::TPC_VIEW_U, inputParameters.m_wirePitchU));
    hitTypeToThickness.insert(HitTypeToFloatMap::value_type(pandora::TPC_VIEW_V, inputParameters.m_wirePitchV));
    hitTypeToThickness.insert(HitTypeToFloatMap::value_type(pandora::TPC_VIEW_W, inputParameters.m_wirePitchW));

    for (const ProtoHit &protoHit : protoHitVector)
    {
        // Mainly dummy parameters
        PandoraApi::CaloHit::Parameters parameters;
        parameters.m_positionVector = pandora::CartesianVector(protoHit.m_x, 0.f, protoHit.m_z);
        parameters.m_expectedDirection = pandora::CartesianVector(0.f, 0.f, 1.f); 
        parameters.m_cellNormalVector = pandora::CartesianVector(0.f, 0.f, 1.f);
        parameters.m_cellGeometry = pandora::RECTANGULAR;
        parameters.m_cellSize0 = 0.5f;
        parameters.m_cellSize1 = inputParameters.m_hitWidth;
        parameters.m_cellThickness = hitTypeToThickness.at(protoHit.m_hitType);
        parameters.m_nCellRadiationLengths = 1.f;
        parameters.m_nCellInteractionLengths = 1.f;
        parameters.m_time = 0.f;
        parameters.m_inputEnergy = protoHit.m_energy;
        parameters.m_mipEquivalentEnergy = 1.f;
        parameters.m_electromagneticEnergy = protoHit.m_energy;
        parameters.m_hadronicEnergy = protoHit.m_energy;
        parameters.m_isDigital = false;
        parameters.m_hitType = protoHit.m_hitType;
        parameters.m_hitRegion = pandora::SINGLE_REGION;
        parameters.m_layer = 0;
        parameters.m_isInOuterSamplingLayer = false;
        parameters.m_pParentAddress = nullptr;

        try
        {
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*pPrimaryPandora, parameters));
        }
        catch (...)
        {
            std::cout << "Unable to make hits" << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------ 

void LoadCell(const Parameters &inputParameters, TiXmlElement *pTiXmlElement, ProtoHitVector &protoHitVectorU, ProtoHitVector &protoHitVectorV, ProtoHitVector &protoHitVectorW)
{
    // ATTN : Geant4 mm, Pandora cm
    const float x(std::atof(pTiXmlElement->Attribute("X"))/10.f);
    const float y(std::atof(pTiXmlElement->Attribute("Y"))/10.f);
    const float z(std::atof(pTiXmlElement->Attribute("Z"))/10.f);
    const pandora::CartesianVector localPosition(x,y,z);

    const float u(YZtoU(localPosition.GetY(), localPosition.GetZ(), inputParameters));
    const float v(YZtoV(localPosition.GetY(), localPosition.GetZ(), inputParameters));
    const float w(localPosition.GetZ());
    const float energy(std::atof(pTiXmlElement->Attribute("Energy")));

    ProtoHit protoHitU, protoHitV, protoHitW;

    protoHitU.m_x = x;
    protoHitU.m_z = u;
    protoHitU.m_energy = energy;
    protoHitU.m_hitType = pandora::TPC_VIEW_U;

    protoHitV.m_x = x;
    protoHitV.m_z = v;
    protoHitV.m_energy = energy;
    protoHitV.m_hitType = pandora::TPC_VIEW_V;

    protoHitW.m_x = x;
    protoHitW.m_z = w;
    protoHitW.m_energy = energy;
    protoHitW.m_hitType = pandora::TPC_VIEW_W;

    protoHitVectorU.push_back(protoHitU);
    protoHitVectorV.push_back(protoHitV);
    protoHitVectorW.push_back(protoHitW);

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LoadMCParticle(TiXmlElement *pTiXmlElement, const pandora::Pandora *const pPrimaryPandora)
{
    PandoraApi::MCParticle::Parameters parameters;

    try
    {
        const float energy(std::atof(pTiXmlElement->Attribute("Energy")));
        const CartesianVector momentum(std::atof(pTiXmlElement->Attribute("MomentumX")), std::atof(pTiXmlElement->Attribute("MomentumY")), std::atof(pTiXmlElement->Attribute("MomentumZ")));
        const CartesianVector vertex(std::atof(pTiXmlElement->Attribute("StartX"))/10.f, std::atof(pTiXmlElement->Attribute("StartY"))/10.f, std::atof(pTiXmlElement->Attribute("StartZ"))/10.f);
        const CartesianVector endpoint(std::atof(pTiXmlElement->Attribute("EndX"))/10.f, std::atof(pTiXmlElement->Attribute("EndY"))/10.f, std::atof(pTiXmlElement->Attribute("EndZ"))/10.f);
        const int particleId(std::atoi(pTiXmlElement->Attribute("PDG")));
        const MCParticleType mcParticleType(pandora::MC_3D);
        const void *pParentAddress(nullptr);

        parameters.m_energy = energy;
        parameters.m_momentum = momentum;
        parameters.m_vertex = vertex;
        parameters.m_endpoint = endpoint;
        parameters.m_particleId = particleId;
        parameters.m_mcParticleType = mcParticleType;
        parameters.m_pParentAddress = pParentAddress;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*pPrimaryPandora, parameters));
    }
    catch (StatusCodeException &statusCodeException)
    {
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float YZtoU(const float y, const float z, const Parameters &parameters)
{
    return (z * std::cos(parameters.m_wireAngleU) - y * std::sin(parameters.m_wireAngleU));
}

//------------------------------------------------------------------------------------------------------------------------------------------ 

float YZtoV(const float y, const float z, const Parameters &parameters)
{
    return (z * std::cos(parameters.m_wireAngleV) - y * std::sin(parameters.m_wireAngleV));
}

//------------------------------------------------------------------------------------------------------------------------------------------ 

void DownsampleHits(const Parameters &inputParameters, ProtoHitVector &protoHitVector)
{
    bool isU(protoHitVector.front().m_hitType == TPC_VIEW_U);
    bool isV(protoHitVector.front().m_hitType == TPC_VIEW_V);
    bool isW(protoHitVector.front().m_hitType == TPC_VIEW_W);

    if (!isU && !isV && !isW)
        throw StopProcessingException("Unexpected hit type");

    const float hitPitch(isU ? inputParameters.m_wirePitchU : isV ? inputParameters.m_wirePitchV : inputParameters.m_wirePitchW);

    if (hitPitch < std::numeric_limits<float>::epsilon())
        throw StopProcessingException("Unfeasibly pitch requested");

    // ATTN : Begin by ordering wire number
    for (ProtoHit &protoHit : protoHitVector)
    {
        if ((isU && protoHit.m_hitType != TPC_VIEW_U) || (isV && protoHit.m_hitType != TPC_VIEW_V) || (isW && protoHit.m_hitType != TPC_VIEW_W))
            throw StopProcessingException("Multiple hit types");

        protoHit.m_z = std::floor((protoHit.m_z + 0.5f * hitPitch) / hitPitch) * hitPitch;
    }

    ProtoHit protoHit1, protoHit2;
    std::sort(protoHitVector.begin(), protoHitVector.end(), SortProtoHits);

    while (IdentifyMerge(inputParameters, protoHitVector, protoHit1, protoHit2))
    {
        ProtoHit mergedHit;

        // ATTN : Merged hit on same wire
        mergedHit.m_z = protoHit1.m_z;
        // ATTN : Energy weighted mean drift position
        mergedHit.m_x = (protoHit1.m_x * protoHit1.m_energy + protoHit2.m_x * protoHit2.m_energy)/(protoHit1.m_energy + protoHit2.m_energy);
        mergedHit.m_energy = protoHit1.m_energy + protoHit2.m_energy;
        mergedHit.m_hitType = protoHit1.m_hitType;

        protoHitVector.erase(std::remove_if(protoHitVector.begin(), protoHitVector.end(), [](const ProtoHit &protoHit) -> bool{return protoHit.m_deleteHit;}), protoHitVector.end());
        protoHitVector.push_back(mergedHit);
        std::sort(protoHitVector.begin(), protoHitVector.end(), SortProtoHits);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------ 

bool IdentifyMerge(const Parameters &inputParameters, ProtoHitVector &protoHitVector, ProtoHit &protoHitA, ProtoHit &protoHitB)
{
    for (int i = 0; i < protoHitVector.size() - 1; i++)
    {
        ProtoHit &protoHit1(protoHitVector.at(i));
        ProtoHit &protoHit2(protoHitVector.at(i+1));

        if (std::fabs(protoHit1.m_z - protoHit2.m_z) < std::numeric_limits<float>::epsilon() && (protoHit1.m_x - protoHit2.m_x < inputParameters.m_hitWidth))
        {
            protoHit1.m_deleteHit = true;
            protoHit2.m_deleteHit = true;
            protoHitA = protoHit1;
            protoHitB = protoHit2;
            return true;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------ 

bool SortProtoHits(const ProtoHit &protoHit1, const ProtoHit &protoHit2)
{
    if (std::fabs(protoHit2.m_z - protoHit1.m_z) > std::numeric_limits<float>::epsilon())
        return protoHit2.m_z > protoHit1.m_z;

    if (std::fabs(protoHit2.m_x - protoHit1.m_x) > std::numeric_limits<float>::epsilon())
        return protoHit2.m_x > protoHit1.m_x;

    return protoHit2.m_energy > protoHit1.m_energy;
}

//------------------------------------------------------------------------------------------------------------------------------------------ 

bool ParseCommandLine(int argc, char *argv[], Parameters &parameters)
{
    if (1 == argc)
        return PrintOptions();

    int c(0);
    std::string recoOption, geometryFileName;

    while ((c = getopt(argc, argv, "r:i:e:n:s:g:pNh")) != -1)
    {
        switch (c)
        {
        case 'r':
            recoOption = optarg;
            break;
        case 'i':
            parameters.m_settingsFile = optarg;
            break;
        case 'e':
            parameters.m_eventFileNameList = optarg;
            break;
        case 'n':
            parameters.m_nEventsToProcess = atoi(optarg);
            break;
        case 's':
            parameters.m_nEventsToSkip = atoi(optarg);
            break;
        case 'g':
            geometryFileName = optarg;
            break;
        case 'p':
            parameters.m_printOverallRecoStatus = true;
            break;
        case 'N':
            parameters.m_shouldDisplayEventNumber = true;
            break;
        case 'h':
        default:
            return PrintOptions();
        }
    }

    TiXmlDocument *pTiXmlDocument = new TiXmlDocument();

    if (!pTiXmlDocument->LoadFile(geometryFileName.c_str()))
        throw StopProcessingException("Invalid geometry xml file");

    const TiXmlHandle xmlDocumentHandle(pTiXmlDocument);
    const TiXmlHandle *pXmlHandle = new TiXmlHandle(xmlDocumentHandle.FirstChildElement().Element());

    try
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(*pXmlHandle, "CenterX", parameters.m_centerX));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(*pXmlHandle, "CenterY", parameters.m_centerY));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(*pXmlHandle, "CenterZ", parameters.m_centerZ));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(*pXmlHandle, "WidthX", parameters.m_widthX));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(*pXmlHandle, "WidthY", parameters.m_widthY));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(*pXmlHandle, "WidthZ", parameters.m_widthZ));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(*pXmlHandle, "WireAngleU", parameters.m_wireAngleU));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(*pXmlHandle, "WireAngleV", parameters.m_wireAngleV));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(*pXmlHandle, "WireAngleW", parameters.m_wireAngleW));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(*pXmlHandle, "WirePitchU", parameters.m_wirePitchU));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(*pXmlHandle, "WirePitchV", parameters.m_wirePitchV));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(*pXmlHandle, "WirePitchW", parameters.m_wirePitchW));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(*pXmlHandle, "HitWidth", parameters.m_hitWidth));
    }
    catch (...)
    {
        throw StopProcessingException("Unable to read geometry component");
    }

    return ProcessRecoOption(recoOption, parameters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LoadXmlElement(const TiXmlElement *pHeadTiXmlElement, float &value, const std::string &name)
{
    const char *attribute(pHeadTiXmlElement->Attribute(name.c_str()));

    if (attribute)
    {
        value = atof(attribute);
        return;
    }

    throw StopProcessingException("Missing geometry parameter " + name);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PrintOptions()
{
    std::cout << std::endl << "./bin/PandoraInterface " << std::endl
              << "    -r RecoOption          (required) [Full, AllHitsCR, AllHitsNu, CRRemHitsSliceCR, CRRemHitsSliceNu, AllHitsSliceCR, AllHitsSliceNu]" << std::endl
              << "    -i Settings            (required) [algorithm description: xml]" << std::endl
              << "    -e EventFileList       (optional) [colon-separated list of files: xml/pndr]" << std::endl
              << "    -n NEventsToProcess    (optional) [no. of events to process]" << std::endl
              << "    -s NEventsToSkip       (optional) [no. of events to skip in first file]" << std::endl
              << "    -u WireAngleU          (optional) [wire angle u, ProtoDUNE assumed]" << std::endl
              << "    -v WireAngleV          (optional) [wire angle v, ProtoDUNE assumed]" << std::endl
              << "    -w WireAngleW          (optional) [wire angle w, ProtoDUNE assumed]" << std::endl
              << "    -p                     (optional) [print status]" << std::endl
              << "    -N                     (optional) [print event numbers]" << std::endl << std::endl;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ProcessRecoOption(const std::string &recoOption, Parameters &parameters)
{
    std::string chosenRecoOption(recoOption);
    std::transform(chosenRecoOption.begin(), chosenRecoOption.end(), chosenRecoOption.begin(), ::tolower);

    if ("full" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = true;
        parameters.m_shouldRunStitching = true;
        parameters.m_shouldRunCosmicHitRemoval = true;
        parameters.m_shouldRunSlicing = true;
        parameters.m_shouldRunNeutrinoRecoOption = true;
        parameters.m_shouldRunCosmicRecoOption = true;
        parameters.m_shouldPerformSliceId = true;
    }
    else if ("allhitscr" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = true;
        parameters.m_shouldRunStitching = true;
        parameters.m_shouldRunCosmicHitRemoval = false;
        parameters.m_shouldRunSlicing = false;
        parameters.m_shouldRunNeutrinoRecoOption = false;
        parameters.m_shouldRunCosmicRecoOption = false;
        parameters.m_shouldPerformSliceId = false;
    }
    else if ("nostitchingcr" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = false;
        parameters.m_shouldRunStitching = false;
        parameters.m_shouldRunCosmicHitRemoval = false;
        parameters.m_shouldRunSlicing = false;
        parameters.m_shouldRunNeutrinoRecoOption = false;
        parameters.m_shouldRunCosmicRecoOption = true;
        parameters.m_shouldPerformSliceId = false;
    } 
    else if ("allhitsnu" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = false;
        parameters.m_shouldRunStitching = false;
        parameters.m_shouldRunCosmicHitRemoval = false;
        parameters.m_shouldRunSlicing = false;
        parameters.m_shouldRunNeutrinoRecoOption = true;
        parameters.m_shouldRunCosmicRecoOption = false;
        parameters.m_shouldPerformSliceId = false;
    }
    else if ("crremhitsslicecr" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = true;
        parameters.m_shouldRunStitching = true;
        parameters.m_shouldRunCosmicHitRemoval = true;
        parameters.m_shouldRunSlicing = true;
        parameters.m_shouldRunNeutrinoRecoOption = false;
        parameters.m_shouldRunCosmicRecoOption = true;
        parameters.m_shouldPerformSliceId = false;
    }
    else if ("crremhitsslicenu" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = true;
        parameters.m_shouldRunStitching = true;
        parameters.m_shouldRunCosmicHitRemoval = true;
        parameters.m_shouldRunSlicing = true;
        parameters.m_shouldRunNeutrinoRecoOption = true;
        parameters.m_shouldRunCosmicRecoOption = false;
        parameters.m_shouldPerformSliceId = false;
    }
    else if ("allhitsslicecr" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = false;
        parameters.m_shouldRunStitching = false;
        parameters.m_shouldRunCosmicHitRemoval = false;
        parameters.m_shouldRunSlicing = true;
        parameters.m_shouldRunNeutrinoRecoOption = false;
        parameters.m_shouldRunCosmicRecoOption = true;
        parameters.m_shouldPerformSliceId = false;
    }
    else if ("allhitsslicenu" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = false;
        parameters.m_shouldRunStitching = false;
        parameters.m_shouldRunCosmicHitRemoval = false;
        parameters.m_shouldRunSlicing = true;
        parameters.m_shouldRunNeutrinoRecoOption = true;
        parameters.m_shouldRunCosmicRecoOption = false;
        parameters.m_shouldPerformSliceId = false;
    }
    else
    {
        std::cout << "LArReco, Unrecognized reconstruction option: " << recoOption << std::endl << std::endl;
        return PrintOptions();
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProcessExternalParameters(const Parameters &parameters, const Pandora *const pPandora)
{
    auto *const pEventReadingParameters = new lar_content::EventReadingAlgorithm::ExternalEventReadingParameters;
    pEventReadingParameters->m_eventFileNameList = parameters.m_eventFileNameList;
    if (parameters.m_nEventsToSkip.IsInitialized()) pEventReadingParameters->m_skipToEvent = parameters.m_nEventsToSkip.Get();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetExternalParameters(*pPandora, "LArEventReading", pEventReadingParameters));

    auto *const pEventSteeringParameters = new lar_content::MasterAlgorithm::ExternalSteeringParameters;
    pEventSteeringParameters->m_shouldRunAllHitsCosmicReco = parameters.m_shouldRunAllHitsCosmicReco;
    pEventSteeringParameters->m_shouldRunStitching = parameters.m_shouldRunStitching;
    pEventSteeringParameters->m_shouldRunCosmicHitRemoval = parameters.m_shouldRunCosmicHitRemoval;
    pEventSteeringParameters->m_shouldRunSlicing = parameters.m_shouldRunSlicing;
    pEventSteeringParameters->m_shouldRunNeutrinoRecoOption = parameters.m_shouldRunNeutrinoRecoOption;
    pEventSteeringParameters->m_shouldRunCosmicRecoOption = parameters.m_shouldRunCosmicRecoOption;
    pEventSteeringParameters->m_shouldPerformSliceId = parameters.m_shouldPerformSliceId;
    pEventSteeringParameters->m_printOverallRecoStatus = parameters.m_printOverallRecoStatus;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetExternalParameters(*pPandora, "LArMaster", pEventSteeringParameters));

#ifdef LIBTORCH_DL
    auto *const pEventSettingsParametersCopy = new lar_content::MasterAlgorithm::ExternalSteeringParameters(*pEventSteeringParameters);
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, pandora::ExternallyConfiguredAlgorithm::SetExternalParameters(*pPandora,
        "LArDLMaster", pEventSettingsParametersCopy));
#endif
}

} // namespace lar_reco
