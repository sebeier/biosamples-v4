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
package uk.ac.ebi.biosamples.certification;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.test.context.junit4.SpringRunner;
import uk.ac.ebi.biosamples.model.certification.Config;
import uk.ac.ebi.biosamples.service.certification.ConfigLoader;

@RunWith(SpringRunner.class)
@SpringBootTest(
    classes = ConfigLoader.class,
    properties = {"job.autorun.enabled=false"})
public class ConfigLoaderTest {
  @Autowired private ConfigLoader configLoader;

  @Test
  public void return_a_valid_config() {
    Config config = configLoader.config;
    assertNotNull(config);
    assertFalse(config.getChecklists().isEmpty());
  }
}
