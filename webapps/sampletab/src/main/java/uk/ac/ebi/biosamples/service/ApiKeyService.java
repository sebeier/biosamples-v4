/*
* Copyright 2019 EMBL - European Bioinformatics Institute
* Licensed under the Apache License, Version 2.0 (the "License"); you may not use this
* file except in compliance with the License. You may obtain a copy of the License at
* http://www.apache.org/licenses/LICENSE-2.0
* Unless required by applicable law or agreed to in writing, software distributed under the
* License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
* CONDITIONS OF ANY KIND, either express or implied. See the License for the
* specific language governing permissions and limitations under the License.
*/
package uk.ac.ebi.biosamples.service;

import java.util.Optional;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.dao.DataAccessException;
import org.springframework.stereotype.Service;
import uk.ac.ebi.biosamples.mongo.model.MongoSampleTabApiKey;
import uk.ac.ebi.biosamples.mongo.repo.MongoSampleTabApiKeyRepository;

@Service
public class ApiKeyService {

  private Logger log = LoggerFactory.getLogger(getClass());

  public static final String BIOSAMPLES = "BioSamples";

  private final MongoSampleTabApiKeyRepository mongoSampleTabApiKeyRepository;

  public ApiKeyService(MongoSampleTabApiKeyRepository mongoSampleTabApiKeyRepository) {
    this.mongoSampleTabApiKeyRepository = mongoSampleTabApiKeyRepository;
  }

  public Optional<String> getDomainForApiKey(String apiKey) throws DataAccessException {
    log.info("getting domain for apikey " + apiKey);
    MongoSampleTabApiKey mongoSampleTabApiKey = mongoSampleTabApiKeyRepository.findOne(apiKey);
    if (mongoSampleTabApiKey == null) return Optional.empty();
    return Optional.ofNullable(mongoSampleTabApiKey.getAapDomain());
  }

  public Optional<String> getUsernameForApiKey(String apiKey) throws DataAccessException {
    log.info("getting domain for apikey " + apiKey);
    MongoSampleTabApiKey mongoSampleTabApiKey = mongoSampleTabApiKeyRepository.findOne(apiKey);
    if (mongoSampleTabApiKey == null) return Optional.empty();
    return Optional.ofNullable(mongoSampleTabApiKey.getUserName());
  }

  public boolean canKeyOwnerEditSource(String keyOwner, String source) {
    if (keyOwner == null || keyOwner.trim().length() == 0) {
      throw new IllegalArgumentException("keyOwner must a sensible string");
    }
    if (source == null || source.trim().length() == 0) {
      throw new IllegalArgumentException("source must be a sensible string");
    }

    if ("BioSamples".toLowerCase().equals(keyOwner.toLowerCase())) {
      // BioSamples key can edit anything
      return true;
    } else if (source.toLowerCase().equals(keyOwner.toLowerCase())) {
      // source key can edit their own samples
      return true;
    } else {
      // deny everyone else
      log.info("Keyowner " + keyOwner + " attempted to access " + source);
      return false;
    }
  }
}
