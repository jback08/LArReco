/**
 *  @file   PandoraSDK/include/MyTemplates/CaloHitWritingAlgorithm.h
 * 
 *  @brief  Header file for the template algorithm class.
 * 
 *  $Log: $
 */
#ifndef CALO_HIT_WRITING_ALGORITHM_H
#define CALO_HIT_WRITING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_reco
{

/**
 *  @brief  CaloHitWritingAlgorithm class
 */
class CaloHitWritingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    // Member variables here
    pandora::StringVector     m_caloHitListNames;    ///< Name of input calo hit list
    std::string               m_trainingOutputFile;  ///< Output file name for training examples
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CaloHitWritingAlgorithm::Factory::CreateAlgorithm() const
{
    return new CaloHitWritingAlgorithm();
}

} // namespace lar_reco

#endif // #ifndef CALO_HIT_WRITING_ALGORITHM_H
