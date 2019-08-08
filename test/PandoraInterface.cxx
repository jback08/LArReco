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
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"
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

    pandora::StringVector eventFileNameVector;
    XmlHelper::TokenizeString(parameters.m_eventFileNameList, eventFileNameVector, ":");

    for (const std::string fileName : eventFileNameVector)
    {
        std::cout << std::endl << "PROCESSING FILE: " << fileName << std::endl;
        // Load Input File
        TiXmlDocument *pTiXmlDocument = new TiXmlDocument();

        if (!pTiXmlDocument->LoadFile(fileName.c_str()))
            std::cerr << pTiXmlDocument->ErrorDesc() << std::endl;

        TiXmlElement *pRunTiXmlElement = pTiXmlDocument->FirstChildElement();
        TiXmlElement *pEventTiXmlElement = pRunTiXmlElement->FirstChildElement();

        while ((nEvents++ < parameters.m_nEventsToProcess) || (0 > parameters.m_nEventsToProcess))
        {
            if (pEventTiXmlElement == 0)
            {
                pTiXmlDocument->Clear();
                delete pTiXmlDocument;
                break;
            }

            if (parameters.m_shouldDisplayEventNumber)
                std::cout << std::endl << "   PROCESSING EVENT: " << (nEvents - 1) << std::endl << std::endl;

            LoadEvent(parameters, pPrimaryPandora, pEventTiXmlElement);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*pPrimaryPandora));
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*pPrimaryPandora));

            pEventTiXmlElement = pEventTiXmlElement->NextSiblingElement();
        }
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
    IntIntMap trackParentId;

    for (TiXmlElement *pSubTiXmlElement = pTiXmlElement->FirstChildElement(); pSubTiXmlElement != NULL; pSubTiXmlElement = pSubTiXmlElement->NextSiblingElement())
    {
        const std::string componentName(pSubTiXmlElement->ValueStr());

        if (componentName == "Cell")
        {
            LoadCell(inputParameters, pSubTiXmlElement, protoHitVectorU, protoHitVectorV, protoHitVectorW);
        }
        else if (componentName == "MCParticle")
        {
            LoadMCParticle(pSubTiXmlElement, pPrimaryPandora, trackParentId);
        }
    }

    // Load MCParticle parent daughter relationships
    for (const auto iter : trackParentId)
    {
        const int particleId(iter.first);
        const int parentId(iter.second);

        if (trackParentId.find(parentId) != trackParentId.end())
        {
            try
            {
                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(*pPrimaryPandora,
                    (void*)((intptr_t)parentId), (void*)((intptr_t)particleId)));
            }
            catch (const pandora::StatusCodeException &)
            {
                std::cout << "LoadEvent - Unable to create mc particle relationship, invalid information supplied " << std::endl;
            }
        }
    }

    if (!protoHitVectorU.empty())
        DownsampleHits(inputParameters, protoHitVectorU);

    if (!protoHitVectorV.empty())
        DownsampleHits(inputParameters, protoHitVectorV);

    if (!protoHitVectorW.empty())
        DownsampleHits(inputParameters, protoHitVectorW);

    ProtoHitVector protoHitVector;
    protoHitVector.insert(protoHitVector.end(), protoHitVectorU.begin(), protoHitVectorU.end());
    protoHitVector.insert(protoHitVector.end(), protoHitVectorW.begin(), protoHitVectorW.end());

    if (!inputParameters.m_dualPhaseMode)
        protoHitVector.insert(protoHitVector.end(), protoHitVectorV.begin(), protoHitVectorV.end());

    HitTypeToFloatMap hitTypeToThickness;
    hitTypeToThickness.insert(HitTypeToFloatMap::value_type(pandora::TPC_VIEW_U, inputParameters.m_wirePitchU));
    hitTypeToThickness.insert(HitTypeToFloatMap::value_type(pandora::TPC_VIEW_V, inputParameters.m_wirePitchV));
    hitTypeToThickness.insert(HitTypeToFloatMap::value_type(pandora::TPC_VIEW_W, inputParameters.m_wirePitchW));

    for (const ProtoHit *pProtoHit : protoHitVector)
    {
        if (pProtoHit->m_energy < inputParameters.m_hitEnergyThreshold)
            continue;

        // Mainly dummy parameters
        lar_content::LArCaloHitParameters parameters;
        parameters.m_positionVector = pandora::CartesianVector(pProtoHit->m_x, 0.f, pProtoHit->m_z);
        parameters.m_expectedDirection = pandora::CartesianVector(0.f, 0.f, 1.f); 
        parameters.m_cellNormalVector = pandora::CartesianVector(0.f, 0.f, 1.f);
        parameters.m_cellGeometry = pandora::RECTANGULAR;
        parameters.m_cellSize0 = 0.5f;
        parameters.m_cellSize1 = inputParameters.m_hitWidth;
        parameters.m_cellThickness = hitTypeToThickness.at(pProtoHit->m_hitType);
        parameters.m_nCellRadiationLengths = 1.f;
        parameters.m_nCellInteractionLengths = 1.f;
        parameters.m_time = 0.f;
        parameters.m_inputEnergy = pProtoHit->m_energy;
        parameters.m_mipEquivalentEnergy = 1.f;
        parameters.m_electromagneticEnergy = pProtoHit->m_energy;
        parameters.m_hadronicEnergy = pProtoHit->m_energy;
        parameters.m_isDigital = false;
        parameters.m_hitType = pProtoHit->m_hitType;
        parameters.m_hitRegion = pandora::SINGLE_REGION;
        parameters.m_layer = 0;
        parameters.m_isInOuterSamplingLayer = false;
        parameters.m_pParentAddress = (void*)((intptr_t)(pProtoHit->m_id));;
        parameters.m_larTPCVolumeId = 0;

        try
        {
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*pPrimaryPandora, parameters));
        }
        catch (...)
        {
            std::cout << "Unable to make hits" << std::endl;
        }

        if (trackParentId.find(pProtoHit->m_mcId) != trackParentId.end())
        {
            try
            {
                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetCaloHitToMCParticleRelationship(*pPrimaryPandora, (void*)((intptr_t)pProtoHit->m_id), (void*)((intptr_t)pProtoHit->m_mcId), 1.f));
            }
            catch (...)
            {
                std::cout << "Unable to make MCParticle - Hit relationships" << std::endl;
            }
        }
        else
        {
            std::cout << "Missing MC particle link to hit" << std::endl;
        }
    }

    for (const ProtoHit *pProtoHit : protoHitVector)
        delete pProtoHit;
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
    const int id(std::atoi(pTiXmlElement->Attribute("Id")));
    const int mcId(std::atoi(pTiXmlElement->Attribute("MCId")));

    ProtoHit *pProtoHitU = new ProtoHit();
    ProtoHit *pProtoHitV = new ProtoHit();
    ProtoHit *pProtoHitW = new ProtoHit();

    pProtoHitU->m_x = x;
    pProtoHitU->m_z = u;
    pProtoHitU->m_energy = energy;
    pProtoHitU->m_hitType = pandora::TPC_VIEW_U;
    pProtoHitU->m_id = id;
    pProtoHitU->m_mcId = mcId;

    pProtoHitV->m_x = x;
    pProtoHitV->m_z = v;
    pProtoHitV->m_energy = energy;
    pProtoHitV->m_hitType = pandora::TPC_VIEW_V;
    pProtoHitV->m_id = id;
    pProtoHitV->m_mcId = mcId;

    pProtoHitW->m_x = x;
    pProtoHitW->m_z = w;
    pProtoHitW->m_energy = energy;
    pProtoHitW->m_hitType = pandora::TPC_VIEW_W;
    pProtoHitW->m_id = id;
    pProtoHitW->m_mcId = mcId;

    protoHitVectorU.push_back(pProtoHitU);
    protoHitVectorV.push_back(pProtoHitV);
    protoHitVectorW.push_back(pProtoHitW);

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LoadMCParticle(TiXmlElement *pTiXmlElement, const pandora::Pandora *const pPrimaryPandora, IntIntMap &trackParentId)
{
    lar_content::LArMCParticleFactory mcParticleFactory;
    lar_content::LArMCParticleParameters parameters;

    try
    {
        const float energy(std::atof(pTiXmlElement->Attribute("Energy")));
        const CartesianVector momentum(std::atof(pTiXmlElement->Attribute("MomentumX")), std::atof(pTiXmlElement->Attribute("MomentumY")), std::atof(pTiXmlElement->Attribute("MomentumZ")));
        const CartesianVector vertex(std::atof(pTiXmlElement->Attribute("StartX"))/10.f, std::atof(pTiXmlElement->Attribute("StartY"))/10.f, std::atof(pTiXmlElement->Attribute("StartZ"))/10.f);
        const CartesianVector endpoint(std::atof(pTiXmlElement->Attribute("EndX"))/10.f, std::atof(pTiXmlElement->Attribute("EndY"))/10.f, std::atof(pTiXmlElement->Attribute("EndZ"))/10.f);
        const int particlePDG(std::atoi(pTiXmlElement->Attribute("PDG")));
        const int particleId(std::atoi(pTiXmlElement->Attribute("Id")));
        const MCParticleType mcParticleType(pandora::MC_3D);

        parameters.m_nuanceCode = 2001;
        parameters.m_energy = energy;
        parameters.m_momentum = momentum;
        parameters.m_vertex = vertex;
        parameters.m_endpoint = endpoint;
        parameters.m_particleId = particlePDG;
        parameters.m_mcParticleType = mcParticleType;
        parameters.m_pParentAddress = (void*)((intptr_t)particleId);

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*pPrimaryPandora, parameters, mcParticleFactory));

        const int parentId(std::atoi(pTiXmlElement->Attribute("ParentId")));
        trackParentId[particleId] = parentId;
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
    bool isU(protoHitVector.front()->m_hitType == TPC_VIEW_U);
    bool isV(protoHitVector.front()->m_hitType == TPC_VIEW_V);
    bool isW(protoHitVector.front()->m_hitType == TPC_VIEW_W);

    if (!isU && !isV && !isW)
        throw StopProcessingException("Unexpected hit type");

    const float hitPitch(isU ? inputParameters.m_wirePitchU : isV ? inputParameters.m_wirePitchV : inputParameters.m_wirePitchW);

    if (hitPitch < std::numeric_limits<float>::epsilon())
        throw StopProcessingException("Unfeasibly pitch requested");

    typedef std::map<int, ProtoHitVector> IntProtoHitVectorMap;
    IntProtoHitVectorMap intProtoHitVectorMap;

    // ATTN : Begin by ordering wire number
    for (ProtoHit *pProtoHit : protoHitVector)
    {
        if ((isU && pProtoHit->m_hitType != TPC_VIEW_U) || (isV && pProtoHit->m_hitType != TPC_VIEW_V) || (isW && pProtoHit->m_hitType != TPC_VIEW_W))
            throw StopProcessingException("Multiple hit types");

        const int wireId(std::floor((pProtoHit->m_z + 0.5f * hitPitch) / hitPitch));
        pProtoHit->m_z = static_cast<float>(wireId) * hitPitch;

        if (intProtoHitVectorMap.find(wireId) != intProtoHitVectorMap.end())
        {
            intProtoHitVectorMap.at(wireId).push_back(pProtoHit);
        }
        else
        {
            ProtoHitVector activeProtoHitVector = {pProtoHit};
            intProtoHitVectorMap.insert(IntProtoHitVectorMap::value_type(wireId, activeProtoHitVector));
        }
    }

    protoHitVector.clear();

    // Merge along x, but only considering hits on same wire at any given time
    for (auto iter : intProtoHitVectorMap)
    {
        ProtoHitVector activeProtoHitVector(iter.second);
        ProtoHit *pProtoHit1(nullptr);
        ProtoHit *pProtoHit2(nullptr);
        std::sort(activeProtoHitVector.begin(), activeProtoHitVector.end(), SortProtoHits);

        while (IdentifyMerge(inputParameters, activeProtoHitVector, pProtoHit1, pProtoHit2))
        {
            ProtoHit *pMergedHit = new ProtoHit();

            // ATTN : Merged hit on same wire
            pMergedHit->m_z = pProtoHit1->m_z;
            // ATTN : Energy weighted mean drift position
            pMergedHit->m_x = (pProtoHit1->m_x * pProtoHit1->m_energy + pProtoHit2->m_x * pProtoHit2->m_energy)/(pProtoHit1->m_energy + pProtoHit2->m_energy);
            pMergedHit->m_energy = pProtoHit1->m_energy + pProtoHit2->m_energy;
            pMergedHit->m_hitType = pProtoHit1->m_hitType;
            pMergedHit->m_id = (pProtoHit1->m_energy > pProtoHit2->m_energy ? pProtoHit1->m_id : pProtoHit2->m_id);
            pMergedHit->m_mcId = (pProtoHit1->m_energy > pProtoHit2->m_energy ? pProtoHit1->m_mcId : pProtoHit2->m_mcId);

            activeProtoHitVector.erase(std::remove(activeProtoHitVector.begin(), activeProtoHitVector.end(), pProtoHit1));
            activeProtoHitVector.erase(std::remove(activeProtoHitVector.begin(), activeProtoHitVector.end(), pProtoHit2));

            delete pProtoHit1;
            delete pProtoHit2;

            // ATTN: Either 1 hit and just push back, or more than one and insert such that x position ordering is preserved
            if (activeProtoHitVector.size() == 0)
            {
                activeProtoHitVector.push_back(pMergedHit);
            }
            else
            {
                bool addedMergedHit(false);

                for (ProtoHitVector::iterator iter2 = activeProtoHitVector.begin(); iter2 != activeProtoHitVector.end(); iter2++)
                {
                    if ((*iter2)->m_x > pMergedHit->m_x)
                    {
                        addedMergedHit = true;
                        activeProtoHitVector.insert(iter2, pMergedHit);
                        break;
                    }
                }

                if (!addedMergedHit)
                    activeProtoHitVector.push_back(pMergedHit);
            }
        }

        protoHitVector.insert(protoHitVector.end(), activeProtoHitVector.begin(), activeProtoHitVector.end());
        activeProtoHitVector.clear();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------ 

bool IdentifyMerge(const Parameters &inputParameters, ProtoHitVector &protoHitVector, ProtoHit *&pProtoHitA, ProtoHit *&pProtoHitB)
{
    for (int i = 0; i < protoHitVector.size() - 1; i++)
    {
        ProtoHit *pProtoHit1(protoHitVector.at(i));
        ProtoHit *pProtoHit2(protoHitVector.at(i+1));

        if (std::fabs(pProtoHit1->m_z - pProtoHit2->m_z) < std::numeric_limits<float>::epsilon() && (std::fabs(pProtoHit1->m_x - pProtoHit2->m_x) < inputParameters.m_hitWidth))
        {
            pProtoHitA = pProtoHit1;
            pProtoHitB = pProtoHit2;
            return true;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------ 

bool SortProtoHits(const ProtoHit *pProtoHit1, const ProtoHit *pProtoHit2)
{
    if (std::fabs(pProtoHit2->m_z - pProtoHit1->m_z) > std::numeric_limits<float>::epsilon())
        return pProtoHit2->m_z > pProtoHit1->m_z;

    if (std::fabs(pProtoHit2->m_x - pProtoHit1->m_x) > std::numeric_limits<float>::epsilon())
        return pProtoHit2->m_x > pProtoHit1->m_x;

    return pProtoHit2->m_energy > pProtoHit1->m_energy;
}

//------------------------------------------------------------------------------------------------------------------------------------------ 

bool ParseCommandLine(int argc, char *argv[], Parameters &parameters)
{
    if (1 == argc)
        return PrintOptions();

    int c(0);
    std::string recoOption, geometryFileName;

    while ((c = getopt(argc, argv, "r:i:e:n:s:g:T:pNDh")) != -1)
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
        case 'T':
            parameters.m_hitEnergyThreshold = atof(optarg);
            break;
        case 'p':
            parameters.m_printOverallRecoStatus = true;
            break;
        case 'N':
            parameters.m_shouldDisplayEventNumber = true;
            break;
        case 'D':
            parameters.m_dualPhaseMode = true;
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
              << "    -N                     (optional) [print event numbers]" << std::endl
              << "    -D                     (optional) [run dual phase, drop v]" << std::endl
              << "    -T threshold           (optional) [hit energy threshold]" << std::endl << std::endl;

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
