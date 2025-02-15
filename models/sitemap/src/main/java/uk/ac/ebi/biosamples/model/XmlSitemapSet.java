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

import java.util.ArrayList;
import java.util.Collection;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElements;
import javax.xml.bind.annotation.XmlRootElement;

/** @author mrelac */
@XmlAccessorType(value = XmlAccessType.NONE)
@XmlRootElement(name = "sitemapset")
public class XmlSitemapSet {
  @XmlElements({@XmlElement(name = "model", type = XmlSitemap.class)})
  private final Collection<XmlSitemap> xmlSitemaps = new ArrayList();

  public void addSitemap(XmlSitemap xmlSitemap) {
    xmlSitemaps.add(xmlSitemap);
  }

  public Collection<XmlSitemap> getXmlSitemaps() {
    return xmlSitemaps;
  }
}
