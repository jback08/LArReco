/**
 *  @file   PandoraSDK/src/MyTemplates/CaloHitWritingAlgorithm.cc
 *
 *  @brief  Implementation of the calo hit writing algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"

#include "CaloHitWritingAlgorithm.h"

using namespace pandora;
using namespace lar_reco;
using namespace lar_content;

StatusCode CaloHitWritingAlgorithm::Run()
{
    for (const std::string listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));

        const bool isU(pCaloHitList->front()->GetHitType() == TPC_VIEW_U ? true : false);
        const bool isV(pCaloHitList->front()->GetHitType() == TPC_VIEW_V ? true : false);
        const bool isW(pCaloHitList->front()->GetHitType() == TPC_VIEW_W ? true : false);

        if (!isU && !isV && !isW)
            return STATUS_CODE_NOT_ALLOWED;

        std::string trainingOutputFileName(m_trainingOutputFile);
        LArMvaHelper::MvaFeatureVector featureVector;

        if (isU) trainingOutputFileName += "_CaloHitListU.txt";
        else if (isV) trainingOutputFileName += "_CaloHitListV.txt";
        else if (isW) trainingOutputFileName += "_CaloHitListW.txt";

        featureVector.push_back(static_cast<double>(pCaloHitList->size()));

        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            int nuanceCode(3000);
            int pdg(-1);

            try
            {
                const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
                nuanceCode = LArMCParticleHelper::GetNuanceCode(LArMCParticleHelper::GetParentMCParticle(pMCParticle));
                pdg = pMCParticle->GetParticleId();
            }
            catch (...)
            {
                continue;
            }

            featureVector.push_back(static_cast<double>(pCaloHit->GetPositionVector().GetX()));
            featureVector.push_back(static_cast<double>(pCaloHit->GetPositionVector().GetY()));
            featureVector.push_back(static_cast<double>(pCaloHit->GetPositionVector().GetZ()));
            featureVector.push_back(static_cast<double>(pdg));
            featureVector.push_back(static_cast<double>(nuanceCode));
       }

        LArMvaHelper::ProduceTrainingExample(trainingOutputFileName, true, featureVector);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloHitWritingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    // Read settings from xml file here
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "CaloHitListNames", m_caloHitListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrainingOutputFileName", m_trainingOutputFile));

    return STATUS_CODE_SUCCESS;
}
