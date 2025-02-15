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
package uk.ac.ebi.biosamples.model.structured;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;
import com.fasterxml.jackson.databind.PropertyNamingStrategy;
import com.fasterxml.jackson.databind.annotation.JsonDeserialize;
import com.fasterxml.jackson.databind.annotation.JsonNaming;
import java.util.*;

@JsonDeserialize(builder = HistologyEntry.Builder.class)
@JsonNaming(PropertyNamingStrategy.SnakeCaseStrategy.class)
@JsonPropertyOrder({"marker", "measurement", "measurement_units", "partner"})
public class HistologyEntry extends StructuredEntry implements Comparable<HistologyEntry> {
  private final StructuredCell marker;
  private final StructuredCell measurement;
  private final StructuredCell measurementUnits;
  private final StructuredCell partner;

  private HistologyEntry(
      StructuredCell marker,
      StructuredCell measurement,
      StructuredCell measurementUnits,
      StructuredCell partner) {
    this.marker = marker;
    this.measurement = measurement;
    this.measurementUnits = measurementUnits;
    this.partner = partner;
  }

  public StructuredCell getMarker() {
    return marker;
  }

  public StructuredCell getMeasurement() {
    return measurement;
  }

  public StructuredCell getMeasurementUnits() {
    return measurementUnits;
  }

  public StructuredCell getPartner() {
    return partner;
  }

  @Override
  public Map<String, StructuredCell> getDataAsMap() {
    // keys must match DataType.HISTOLOGY.getHeaders()
    Map<String, StructuredCell> keyToStructuredCellMap = new HashMap<>();

    keyToStructuredCellMap.put("Marker", marker);
    keyToStructuredCellMap.put("Measurement", measurement);
    keyToStructuredCellMap.put("Measurement Units", measurementUnits);
    keyToStructuredCellMap.put("Partner", partner);

    return keyToStructuredCellMap;
  }

  public int compareTo(HistologyEntry other) {
    if (other == null) {
      return 1;
    }

    int cmp = this.marker.compareTo(other.marker);
    if (cmp != 0) {
      return cmp;
    }
    cmp = this.measurement.compareTo(other.measurement);
    if (cmp != 0) {
      return cmp;
    }
    cmp = this.measurementUnits.compareTo(other.measurementUnits);
    if (cmp != 0) {
      return cmp;
    }
    return this.partner.compareTo(other.partner);
  }

  @Override
  public boolean equals(Object o) {
    if (this == o) return true;
    if (o instanceof HistologyEntry) {
      HistologyEntry other = (HistologyEntry) o;
      return Objects.equals(this.getMarker(), other.getMarker())
          && Objects.equals(this.getMeasurement(), other.getMeasurement())
          && Objects.equals(this.getMeasurementUnits(), other.getMeasurementUnits())
          && Objects.equals(this.getPartner(), other.getPartner());
    }

    return false;
  }

  @Override
  public int hashCode() {
    return Objects.hash(marker, measurement, measurementUnits, partner);
  }

  @JsonNaming(PropertyNamingStrategy.SnakeCaseStrategy.class)
  public static class Builder {
    private StructuredCell marker;
    private StructuredCell measurement;
    private StructuredCell measurementUnits;
    private StructuredCell partner;

    @JsonCreator
    public Builder() {
      // Start with empty builder
    }

    @JsonProperty
    public Builder withMarker(StructuredCell marker) {
      this.marker = marker;
      return this;
    }

    @JsonProperty
    public Builder withMeasurement(StructuredCell measurement) {
      this.measurement = measurement;
      return this;
    }

    @JsonProperty
    public Builder withMeasurementUnits(StructuredCell measurementUnits) {
      this.measurementUnits = measurementUnits;
      return this;
    }

    @JsonProperty
    public Builder withPartner(StructuredCell partner) {
      this.partner = partner;
      return this;
    }

    public HistologyEntry build() {
      return new HistologyEntry(this.marker, this.measurement, this.measurementUnits, this.partner);
    }
  }
}
