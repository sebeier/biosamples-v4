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
package uk.ac.ebi.biosamples.exceptions;

public class ConflictingSampleTabOwnershipException extends SampleTabException {

  private static final long serialVersionUID = -1504945560846665587L;
  public final String sampleAccession;
  public final String originalSubmission;
  public final String newSubmission;

  public ConflictingSampleTabOwnershipException(
      String sampleAccession, String originalSubmission, String newSubmission) {
    super("Accession " + sampleAccession + " was previously described in " + originalSubmission);
    this.sampleAccession = sampleAccession;
    this.originalSubmission = originalSubmission;
    this.newSubmission = newSubmission;
  }
}
