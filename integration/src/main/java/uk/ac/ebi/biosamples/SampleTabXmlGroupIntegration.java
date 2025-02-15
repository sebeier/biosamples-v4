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
package uk.ac.ebi.biosamples;

import java.net.URI;
import java.nio.charset.StandardCharsets;
import java.util.Optional;
import java.util.Scanner;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.boot.web.client.RestTemplateBuilder;
import org.springframework.core.annotation.Order;
import org.springframework.hateoas.Resource;
import org.springframework.http.MediaType;
import org.springframework.http.RequestEntity;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Component;
import org.springframework.web.client.RestOperations;
import org.springframework.web.util.UriComponentsBuilder;
import uk.ac.ebi.biosamples.client.BioSamplesClient;
import uk.ac.ebi.biosamples.model.Sample;

@Component
@Order(6)
public class SampleTabXmlGroupIntegration extends AbstractIntegration {

  private Logger log = LoggerFactory.getLogger(this.getClass());

  private final RestOperations restTemplate;

  private final URI uri;

  public SampleTabXmlGroupIntegration(
      RestTemplateBuilder restTemplateBuilder,
      IntegrationProperties integrationProperties,
      BioSamplesClient client) {
    super(client);
    this.restTemplate = restTemplateBuilder.build();

    uri =
        UriComponentsBuilder.fromUri(integrationProperties.getBiosampleSubmissionUriSampleTab())
            .pathSegment("api", "v2", "source", "biosamples", "group")
            .queryParam("apikey", integrationProperties.getLegacyApiKey())
            .build()
            .toUri();
  }

  @Override
  protected void phaseOne() {

    runCallableOnResource(
        "/SAMEG123_unaccession.xml",
        sampleTabString -> {
          log.info("POSTing to " + uri);
          RequestEntity<String> request =
              RequestEntity.post(uri)
                  .contentType(MediaType.APPLICATION_XML)
                  .accept(MediaType.TEXT_PLAIN)
                  .body(sampleTabString);
          ResponseEntity<String> response = restTemplate.exchange(request, String.class);
          String accession = response.getBody();
          // check at the right URLs with GET to make sure all
          // arrived
          Optional<Resource<Sample>> group = client.fetchSampleResource(accession);
          if (!group.isPresent()) {
            throw new RuntimeException("Unable to retrieve group " + accession);
          }
          if (group.get().getContent().getAttributes().size() == 0) {
            throw new RuntimeException("No attributes on group " + accession);
          }
          if (group.get().getContent().getRelationships().size() == 0) {
            throw new RuntimeException("No relationships on group " + accession);
          }
        });

    URI putUri = UriComponentsBuilder.fromUri(uri).pathSegment("SAMEG123").build().toUri();
    runCallableOnResource(
        "/SAMEG123.xml",
        sampleTabString -> {
          log.info("PUTing to " + putUri);
          RequestEntity<String> request =
              RequestEntity.put(putUri)
                  .contentType(MediaType.APPLICATION_XML)
                  .accept(MediaType.TEXT_PLAIN)
                  .body(sampleTabString);
          ResponseEntity<String> response = restTemplate.exchange(request, String.class);
          String accession = response.getBody();
          // check at the right URLs with GET to make sure all
          // arrived
          Optional<Resource<Sample>> group = client.fetchSampleResource(accession);
          if (!group.isPresent()) {
            throw new RuntimeException("Unable to retrieve group " + accession);
          }
          if (group.get().getContent().getAttributes().size() == 0) {
            throw new RuntimeException("No attributes on group " + accession);
          }
          if (group.get().getContent().getRelationships().size() == 0) {
            throw new RuntimeException("No relationships on group " + accession);
          }
        });
  }

  @Override
  protected void phaseTwo() {}

  @Override
  protected void phaseThree() {}

  @Override
  protected void phaseFour() {
    // TODO Auto-generated method stub

  }

  @Override
  protected void phaseFive() {
    // TODO Auto-generated method stub

  }

  private interface Callback {
    void callback(String sampleTabString);
  }

  private void runCallableOnResource(String resource, Callback callback) {

    Scanner scanner = null;
    String xmlString;

    try {
      scanner =
          new Scanner(
              this.getClass().getResourceAsStream(resource), StandardCharsets.UTF_8.toString());
      xmlString = scanner.useDelimiter("\\A").next();
    } finally {
      if (scanner != null) {
        scanner.close();
      }
    }

    log.trace("sending legacy xml submission \n" + xmlString);

    if (xmlString != null) {
      callback.callback(xmlString);
    }
  }
}
