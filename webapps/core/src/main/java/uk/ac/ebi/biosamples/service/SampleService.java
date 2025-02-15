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

import java.time.Instant;
import java.util.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.http.HttpStatus;
import org.springframework.stereotype.Service;
import org.springframework.web.bind.annotation.ResponseStatus;
import uk.ac.ebi.biosamples.model.Autocomplete;
import uk.ac.ebi.biosamples.model.Sample;
import uk.ac.ebi.biosamples.model.StaticViewWrapper;
import uk.ac.ebi.biosamples.model.filter.Filter;
import uk.ac.ebi.biosamples.mongo.model.MongoRelationship;
import uk.ac.ebi.biosamples.mongo.model.MongoSample;
import uk.ac.ebi.biosamples.mongo.repo.MongoSampleRepository;
import uk.ac.ebi.biosamples.mongo.service.MongoAccessionService;
import uk.ac.ebi.biosamples.mongo.service.MongoSampleToSampleConverter;
import uk.ac.ebi.biosamples.mongo.service.SampleToMongoSampleConverter;
import uk.ac.ebi.biosamples.mongo.service.SampleToMongoSampleStructuredDataCentricConverter;
import uk.ac.ebi.biosamples.solr.service.SolrSampleService;

/**
 * Service layer business logic for centralising repository access and conversions between different
 * controller. Use this instead of linking to repositories directly.
 *
 * @author faulcon
 */
@Service
public class SampleService {
  private static final String VALIDATION_MESSAGE =
      "Only Sample name, sample accession and sample structured data can be provided through this API";
  private static final String NO_STRUCTURED_DATA_IS_PROVIDED = "No structured data is provided";
  private static final String NCBI_IMPORT_DOMAIN = "self.BiosampleImportNCBI";
  private static final String ENA_IMPORT_DOMAIN = "self.BiosampleImportENA";
  private static Logger log = LoggerFactory.getLogger(SampleService.class);

  // TODO use constructor injection
  @Qualifier("SampleAccessionService")
  @Autowired
  private MongoAccessionService mongoAccessionService;

  @Qualifier("GroupAccessionService")
  @Autowired
  private MongoAccessionService groupAccessionService;

  @Autowired private MongoSampleRepository mongoSampleRepository;
  @Autowired private MongoSampleToSampleConverter mongoSampleToSampleConverter;
  @Autowired private SampleToMongoSampleConverter sampleToMongoSampleConverter;
  @Autowired private SampleToMongoSampleStructuredDataCentricConverter structuredDataConverter;
  @Autowired private SampleValidator sampleValidator;
  @Autowired private SolrSampleService solrSampleService;
  @Autowired private SampleReadService sampleReadService;
  @Autowired private MessagingService messagingSerivce;

  /**
   * Throws an IllegalArgumentException of no sample with that accession exists
   *
   * @param accession the sample accession
   * @return
   * @throws IllegalArgumentException
   */
  public Optional<Sample> fetch(
      String accession, Optional<List<String>> curationDomains, String curationRepo) {
    StaticViewWrapper.StaticView staticView =
        StaticViewWrapper.getStaticView(curationDomains.orElse(null), curationRepo);
    return sampleReadService.fetch(accession, curationDomains, staticView);
  }

  public Autocomplete getAutocomplete(
      String autocompletePrefix, Collection<Filter> filters, int noSuggestions) {
    return solrSampleService.getAutocomplete(autocompletePrefix, filters, noSuggestions);
  }

  public boolean beforeStore(Sample sample) {
    return beforeStoreCheck(sample);
  }

  private boolean beforeStoreCheck(Sample sample) {
    boolean firstTimeMetadataAdded;

    final String domain = sample.getDomain();

    if (isPipelineEnaOrNcbiDomain(domain))
      firstTimeMetadataAdded = false; // imported sample - never submitted first time to BSD
    else {
      firstTimeMetadataAdded = isFirstTimeMetadataAddedForNonImportedSamples(sample);
    }

    if (firstTimeMetadataAdded) log.trace("First time metadata added");

    return firstTimeMetadataAdded;
  }

  private boolean isFirstTimeMetadataAddedForNonImportedSamples(Sample sample) {
    boolean firstTimeMetadataAdded = true;

    if (sample.hasAccession()) {
      MongoSample mongoOldSample = mongoSampleRepository.findOne(sample.getAccession());
      if (mongoOldSample != null) {
        firstTimeMetadataAdded = isFirstTimeMetadataAdded(firstTimeMetadataAdded, mongoOldSample);
      }
    }

    return firstTimeMetadataAdded;
  }

  private boolean isPipelineEnaOrNcbiDomain(String domain) {
    return isPipelineEnaDomain(domain) || isPipelineNcbiDomain(domain);
  }

  private boolean isFirstTimeMetadataAdded(
      boolean firstTimeMetadataAdded, MongoSample mongoOldSample) {
    Sample oldSample = mongoSampleToSampleConverter.convert(mongoOldSample);

    if (oldSample.getTaxId() != null && oldSample.getTaxId() > 0) {
      firstTimeMetadataAdded = false;
    }

    if (oldSample.getAttributes().size() > 0) {
      firstTimeMetadataAdded = false;
    }

    if (oldSample.getRelationships().size() > 0) {
      firstTimeMetadataAdded = false;
    }

    if (oldSample.getPublications().size() > 0) {
      firstTimeMetadataAdded = false;
    }

    if (oldSample.getContacts().size() > 0) {
      firstTimeMetadataAdded = false;
    }

    if (oldSample.getOrganizations().size() > 0) {
      firstTimeMetadataAdded = false;
    }

    return firstTimeMetadataAdded;
  }

  // because the fetchUsing caches the sample, if an updated version is stored, we need to make
  // sure
  // that any cached version
  // is removed.
  // Note, pages of samples will not be cache busted, only single-accession sample retrieval
  // @CacheEvict(cacheNames=WebappProperties.fetchUsing, key="#result.accession")
  public Sample store(Sample sample, boolean isFirstTimeMetadataAdded, String authProvider) {
    return store(sample, false, isFirstTimeMetadataAdded, authProvider);
  }

  public Sample store(
      Sample sample, boolean sampleGroup, boolean isFirstTimeMetadataAdded, String authProvider) {
    Collection<String> errors = sampleValidator.validate(sample);
    if (!errors.isEmpty()) {
      log.error("Sample validation failed : {}", errors);
      throw new SampleValidationException(String.join("|", errors));
    }

    if (sample.hasAccession()) {
      MongoSample mongoOldSample = mongoSampleRepository.findOne(sample.getAccession());
      List<String> existingRelationshipTargets = new ArrayList<>();

      if (mongoOldSample != null) {
        Sample oldSample = mongoSampleToSampleConverter.convert(mongoOldSample);
        existingRelationshipTargets =
            getExistingRelationshipTargets(sample.getAccession(), mongoOldSample);

        sample =
            compareWithExistingAndUpdateSample(
                sample, oldSample, isFirstTimeMetadataAdded, authProvider);
      } else {
        log.error("Trying to update sample not in database, accession: {}", sample.getAccession());
      }

      MongoSample mongoSample = sampleToMongoSampleConverter.convert(sample);

      mongoSample = mongoSampleRepository.save(mongoSample);
      sample = mongoSampleToSampleConverter.convert(mongoSample);

      // send a message for storage and further processing, send relationship targets to
      // identify
      // deleted relationships
      messagingSerivce.fetchThenSendMessage(sample.getAccession(), existingRelationshipTargets);
    } else {
      if (sampleGroup) {
        sample = groupAccessionService.generateAccession(sample);
      } else {
        sample = mongoAccessionService.generateAccession(sample);
      }
      messagingSerivce.fetchThenSendMessage(sample.getAccession());
    }

    // do a fetch to return it with accession, curation objects, inverse relationships
    return fetch(sample.getAccession(), Optional.empty(), null).get();
  }

  public Sample storeSampleStructuredData(Sample newSample, String authProvider) {
    validateSampleContentsForStructuredDataPatching(newSample);

    MongoSample mongoOldSample = mongoSampleRepository.findOne(newSample.getAccession());

    if (mongoOldSample != null) {
      newSample = makeNewSample(newSample, mongoSampleToSampleConverter.convert(mongoOldSample));
    } else {
      log.error(
          "Trying to update newSample not in database, accession: {}", newSample.getAccession());
    }

    MongoSample mongoSample = structuredDataConverter.convert(newSample, authProvider);
    mongoSample = mongoSampleRepository.save(mongoSample);
    newSample = mongoSampleToSampleConverter.convert(mongoSample);

    // return the newSample in case we have modified it i.e accessioned
    // do a fetch to return it with curation objects and inverse relationships
    return fetch(newSample.getAccession(), Optional.empty(), null).get();
  }

  private void validateSampleContentsForStructuredDataPatching(Sample newSample) {
    assert newSample.getData() != null;

    final String domain = newSample.getDomain();

    if (!(newSample.getData().size() > 0)) {
      throw new SampleValidationException(NO_STRUCTURED_DATA_IS_PROVIDED);
    }

    if (newSample.getAttributes() != null && newSample.getAttributes().size() > 0) {
      throw new SampleValidationException(VALIDATION_MESSAGE);
    }

    if (newSample.getExternalReferences() != null && newSample.getExternalReferences().size() > 0) {
      throw new SampleValidationException(VALIDATION_MESSAGE);
    }

    if (newSample.getRelationships() != null && newSample.getRelationships().size() > 0) {
      throw new SampleValidationException(VALIDATION_MESSAGE);
    }

    if (newSample.getContacts() != null && newSample.getContacts().size() > 0) {
      throw new SampleValidationException(VALIDATION_MESSAGE);
    }

    if (newSample.getPublications() != null && newSample.getPublications().size() > 0) {
      throw new SampleValidationException(VALIDATION_MESSAGE);
    }

    if (domain != null && domain.length() > 0) {
      throw new SampleValidationException(VALIDATION_MESSAGE);
    }

    if (!newSample.hasAccession()) {
      throw new SampleValidationException("Sample doesn't have an accession");
    }
  }

  private Sample makeNewSample(Sample newSample, Sample oldSample) {
    return Sample.Builder.fromSample(oldSample)
        .withData(newSample.getData())
        .withUpdate(Instant.now())
        .build();
  }

  public boolean searchSampleByDomainAndName(final String domain, final String name) {
    return mongoSampleRepository.findByDomainAndName(domain, name).size() > 0;
  }

  public void validateSample(Map sampleAsMap) {
    List<String> errors = sampleValidator.validate(sampleAsMap);
    StringBuilder sb = new StringBuilder();
    if (errors.size() > 0) {
      for (String error : errors) {
        sb.append(error).append("; ");
      }

      throw new SampleValidationException(sb.toString());
    }
  }

  public boolean isNotExistingAccession(String accession) {
    return mongoSampleRepository.findOne(accession) == null;
  }

  @ResponseStatus(HttpStatus.BAD_REQUEST)
  public static class SampleValidationException extends RuntimeException {
    private static final long serialVersionUID = -7937033504537036300L;

    public SampleValidationException(String message) {
      super(message);
    }
  }

  private List<String> getExistingRelationshipTargets(
      String accession, MongoSample mongoOldSample) {
    List<String> oldRelationshipTargets = new ArrayList<>();
    for (MongoRelationship relationship : mongoOldSample.getRelationships()) {
      if (relationship.getSource().equals(accession)) {
        oldRelationshipTargets.add(relationship.getTarget());
      }
    }

    return oldRelationshipTargets;
  }

  private Sample compareWithExistingAndUpdateSample(
      Sample sampleToUpdate,
      Sample oldSample,
      boolean isFirstTimeMetadataAdded,
      String authProvider) {
    // compare with existing version and check what fields have changed
    if (sampleToUpdate.equals(oldSample)) {
      log.info("New sample is similar to the old sample, accession: {}", oldSample.getAccession());
    }

    // Keep the create date as existing sample -- earlier
    // 13/01/2020 - if the sample has a date, acknowledge it. It can be the actual create date
    // from
    // NCBI, ENA.

    return Sample.Builder.fromSample(sampleToUpdate)
        .withCreate(defineCreateDate(sampleToUpdate, oldSample, authProvider))
        .withSubmitted(
            defineSubmittedDate(sampleToUpdate, oldSample, isFirstTimeMetadataAdded, authProvider))
        // .withUpdate(defineUpdatedDate(sampleToUpdate, oldSample))
        .build();
  }

  private Instant defineUpdatedDate(Sample sampleToUpdate, Sample oldSample) {
    if (isPipelineEnaOrNcbiDomain(sampleToUpdate.getDomain())) {
      return sampleToUpdate.getUpdate() != null
          ? sampleToUpdate.getUpdate()
          : oldSample.getUpdate() != null ? oldSample.getUpdate() : Instant.now();
    } else {
      return Instant.now();
    }
  }

  private Instant defineCreateDate(
      final Sample sampleToUpdate, final Sample oldSample, String authProvider) {
    final String domain = sampleToUpdate.getDomain();

    if (!isWebinAuthorization(authProvider)) {
      if (isPipelineNcbiDomain(domain)) {
        return sampleToUpdate.getCreate() != null
            ? sampleToUpdate.getCreate()
            : (oldSample.getCreate() != null ? oldSample.getCreate() : oldSample.getUpdate());
      } else if (isPipelineEnaDomain(domain)) {
        return (oldSample.getCreate() != null) ? oldSample.getCreate() : sampleToUpdate.getCreate();
      }
    }

    return (oldSample.getCreate() != null ? oldSample.getCreate() : oldSample.getUpdate());
  }

  public boolean isWebinAuthorization(String authProviderIdentifier) {
    return authProviderIdentifier != null && authProviderIdentifier.equalsIgnoreCase("WEBIN");
  }

  private boolean isPipelineEnaDomain(String domain) {
    if (domain == null) return false;
    return domain.equalsIgnoreCase(ENA_IMPORT_DOMAIN);
  }

  private boolean isPipelineNcbiDomain(String domain) {
    if (domain == null) return false;
    return domain.equalsIgnoreCase(NCBI_IMPORT_DOMAIN);
  }

  private Instant defineSubmittedDate(
      final Sample sampleToUpdate,
      final Sample oldSample,
      boolean isFirstTimeMetadataAdded,
      String authProvider) {
    final String domain = sampleToUpdate.getDomain();

    if (isPipelineNcbiDomain(domain)) {
      return sampleToUpdate.getSubmitted() != null
          ? sampleToUpdate.getSubmitted()
          : (oldSample.getSubmitted() != null ? oldSample.getSubmitted() : oldSample.getCreate());
    } else if (isPipelineEnaDomain(domain)) {
      return (oldSample.getSubmitted() != null)
          ? oldSample.getSubmitted()
          : sampleToUpdate.getSubmitted();
    } else {
      if (isFirstTimeMetadataAdded) {
        return sampleToUpdate.getSubmitted();
      } else {
        return oldSample.getSubmitted() != null
            ? oldSample.getSubmitted()
            : (oldSample.getCreate() != null ? oldSample.getCreate() : oldSample.getUpdate());
      }
    }
  }

  /*
  //this code recursively follows relationships
  //TODO finish
  public SortedSet<Sample> getRelated(Sample sample, String relationshipType) {
  	Queue<String> toCheck = new LinkedList<>();
  	Set<String> checked = new HashSet<>();
  	Collection<Sample> related = new TreeSet<>();
  	toCheck.add(sample.getAccession());
  	while (!toCheck.isEmpty()) {
  		String accessionToCheck = toCheck.poll();
  		checked.add(accessionToCheck);
  		Sample sampleToCheck = sampleReadService.fetchUsing(accessionToCheck);
  		related.add(sampleToCheck);
  		for (Relationship rel : sampleToCheck.getRelationships()) {
  			if (relationshipType == null || relationshipType.equals(rel.getType())) {
  				if (!checked.contains(rel.getSource()) && toCheck.contains(rel.getSource())) {
  					toCheck.add(rel.getSource());
  				}
  				if (!checked.contains(rel.getTarget()) && toCheck.contains(rel.getTarget())) {
  					toCheck.add(rel.getTarget());
  				}
  			}
  		}
  	}
  	related.remove(sample);
  	return related;
  }
  */
}
