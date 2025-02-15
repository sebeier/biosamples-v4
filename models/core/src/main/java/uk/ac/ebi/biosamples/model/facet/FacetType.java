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
package uk.ac.ebi.biosamples.model.facet;

import com.fasterxml.jackson.annotation.JsonValue;
import java.lang.reflect.InvocationTargetException;

public enum FacetType {
  ATTRIBUTE_FACET("attribute", AttributeFacet.Builder.class),
  DATE_RANGE_FACET("date range", DateRangeFacet.Builder.class),
  INVERSE_RELATION_FACET("inverse relation", InverseRelationFacet.Builder.class),
  RELATION_FACET("relation", RelationFacet.Builder.class),
  EXTERNAL_REFERENCE_DATA_FACET(
      "external reference data", ExternalReferenceDataFacet.Builder.class),
  NO_TYPE("", null);

  private String name;
  private Class<? extends Facet.Builder> associatedClass;

  FacetType(String name, Class<? extends Facet.Builder> associatedClass) {
    this.name = name;
    this.associatedClass = associatedClass;
  }

  public Facet.Builder getBuilderForLabelAndCount(String facetLabel, Long facetCount) {
    try {
      return this.associatedClass
          .getConstructor(String.class, Long.class)
          .newInstance(facetLabel, facetCount);
    } catch (NoSuchMethodException
        | IllegalAccessException
        | InstantiationException
        | InvocationTargetException e) {
      throw new RuntimeException("FacetType " + this + " does not provide a proper builder");
    }
  }

  @JsonValue
  public String getFacetName() {
    return this.name.toLowerCase();
  }
}
