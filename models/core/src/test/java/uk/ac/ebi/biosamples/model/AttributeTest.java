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
package uk.ac.ebi.biosamples.model;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.assertThat;

import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;
import java.nio.charset.StandardCharsets;
import org.junit.Test;
import uk.ac.ebi.biosamples.service.HttpOlsUrlResolutionService;

// @RunWith(SpringRunner.class)
// @JsonTest
public class AttributeTest {

  @Test
  public void test_getIriOls_method_returns_null_if_invalid_iri_is_provided() {
    Attribute testAttribute =
        Attribute.build("WrongIRIAttributeKey", "WrongIRIAttributeValue", "Something else", null);

    HttpOlsUrlResolutionService httpOlsUrlResolutionService = new HttpOlsUrlResolutionService();

    assertThat(httpOlsUrlResolutionService.getIriOls(testAttribute.getIri()), nullValue());
  }

  @Test
  public void test_getIriOls_method_returns_iri_if_valid_iri_is_provided()
      throws UnsupportedEncodingException {
    String iri = "http://purl.obolibrary.org/obo/NCBITaxon_291302";
    Attribute testAttribute = Attribute.build("Organism", "Miniopterus natalensis", iri, null);

    HttpOlsUrlResolutionService httpOlsUrlResolutionService = new HttpOlsUrlResolutionService();

    System.out.println(httpOlsUrlResolutionService.getIriOls(testAttribute.getIri()));
    assertThat(
        httpOlsUrlResolutionService.getIriOls(testAttribute.getIri()),
        allOf(
            endsWith("NCBITaxon_291302"),
            startsWith("http://www.ebi.ac.uk/ols/terms?iri="),
            containsString(URLEncoder.encode(iri, StandardCharsets.UTF_8.toString()))));
  }

  @Test
  public void test_getIriOls_method_returns_null_for_a_query() {
    String curie = "CL:0000451";
    Attribute testAttribute = Attribute.build("Cell type", "dendritic cell", curie, null);

    HttpOlsUrlResolutionService httpOlsUrlResolutionService = new HttpOlsUrlResolutionService();

    assertThat(httpOlsUrlResolutionService.getIriOls(testAttribute.getIri()), nullValue());
  }
}
