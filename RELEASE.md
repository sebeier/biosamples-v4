# 4.0.0

Version v4.0.0 represents a re-architecture and re-engineering of the BioSamples software stack. It is now based on the Java [Spring-Boot](https://projects.spring.io/spring-boot) framework, utilising [MongoDB](https://www.mongodb.com) for storage and [Solr](https://lucene.apache.org/solr) for indexing and search. It tries to follow up-to-date web standards and conventions, while remaining backwards compatible. This will also give us a strong and stable foundation to build more features and improvements from, more reliably and more rapidly.

Highlights include:

* Submissions and updates will be available immediately via accession, and will be available via search within a few minutes or less. There is also improved handling of submissions and updates, with fewer errors and better feedback about any problems.
* Integration with [EBI AAP](https://aap.tsi.ebi.ac.uk) for login management and access to pre-publication samples, including use of [ELIXIR AAI](https://www.elixir-europe.org/services/compute/aai) single sign-on accounts.
* Separation of submitted sample information from curation of that information, including the ability for 3rd party (re-)curation of samples. Please contact us if you would be interested in more information and/or to supply curation information.
* Improved handling of non-alphanumeric characters in attribute types e.g. "geographic location (country and/or sea)"
* Improved faceting allowing selection of multiple values within same facet, fixed re-use and re-distribution of search URLs. This will be expanded in future with additional facet types where appropriate.
* Support and recommend the use of [content negotiation](https://developer.mozilla.org/en-US/docs/Web/HTTP/Content_negotiation) to accessing multiple formats at the same URIs. In addition to the content (HTML vs XML vs JSON) this also supports [compression](https://developer.mozilla.org/en-US/docs/Web/HTTP/Compression) and [caching](https://developer.mozilla.org/en-US/docs/Web/HTTP/Caching) through standard mechanisms.
* Java client using Spring, and a Spring-Boot starter module for easy use. This is used by BioSamples internally and other teams at EMBL-EBI, so is high performance and battle tested.
* Containerisation using Docker and Docker-Compose, which makes it easier to run a local version for client development or for local storage of sample information.

## Data content

* Ontology terms Numeric tax IDs (e.g. 9606) and short ontology terms (e.g. PATO:0000384) are being replaced with full IRIs (e.g. <http://purl.obolibrary.org/obo/NCBITaxon_9606> and <http://purl.obolibrary.org/obo/PATO_0000384> ) in many places, eventually everywhere.
* Groups will continue to exist for backwards compatibility purposes. However, we are investigating future development to reduce or remove many of these in favour of alternatives such as filtering samples by external link, or delegating grouping of samples to other EMBL-EBI archives such as [BioStudies](https://www.ebi.ac.uk/biostudies).

## JSON `/biosamples`

This is the preferred API for use, and uses the same URIs as the HTML pages, and utilising content negotiation to provide a JSON response. This is designed as a [hypermedia as the engine of application state (HATEOS) API](https://en.wikipedia.org/wiki/Hypertext_Application_Language) and therefore we recommend users do not use specific URLs but rather follow relationships between API endpoints, much like a user would use links between HTML pages. It is similar to the `/biosamples/api` JSON format, with a few critical differences:

* added *release* in full ISO 8601 format including time. The backwards-compatible *releaseDate* exists but should be considered deprecated and will be removed in a future release.
* added *update* in full ISO 8601 format including time. The backwards-compatible *updateDate* exists but should be considered deprecated and will be removed in a future release.
* removed *description* as a separate field, is now available as a *characteristic*. 
* remove **relations** rel link; equivalent information is now embedded in sample in *relationships* and *externalReferences* lists.
* remove **sample** rel link; with relations now embedded, this link serves no purpose.
* added **curationLinks** rel link.
* ordering may be different.
* fields are not displayed if empty or null.
* characteristic names accurately reflect what was submitted and may now be multiple words and may include non alphanumeric characters (e.g brackets, greek letters, etc). In the `/biosamples/api` responses characteristic names were always camelCased and with non-alphanumeric characters removed.
* external references directly embedded in the samples and the groups.
  
## XML `/biosamples/xml`

We are maintaining this for backwards compatibility. Later in 2018 we will be consulting about future development of this API, particularly in the context of the improved JSON `/biosamples` API using content negotiation and several long-standing issues with limitations arising from the XML schema in use.

* XML element **TermSourceREF** element **Name** and element **URI** are removed.
* XML element **Property** attributes characteristic and comment always false.
* elements and attributes may be in different order.
* allows only one IRI on attributes, so in rare cases of multiple IRIs will not be complete.
* Query parameter `query` has now a default value of * if none is provided.
* Query parameter `sort` is ignored for the search, due to undefined behaviour and lack of usage.
  
## JSON `/biosamples/api`

This API should be considered **deprecated** and we will aim to remove it by 2019. Any users of this should move to using the `/biosamples` URIs to retrieve JSON representations with an improved schema via content negotiation. Further announcements will be made in future for specific updates and deadlines.

* ordering may be different from previous versions, and is not guaranteed for future versions.
* fields are not displayed if empty or null.
* `/api/externallinksrelations/{id}/sample` and `/api/externallinksrelations/{id}/group` are removed due to lack of usage.
* fixed *externalReferences* and *publications* to be nested objects and not JSON strings.