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
package uk.ac.ebi.biosamples.ebeye.base;

import com.mongodb.*;
import java.io.File;
import java.util.*;
import java.util.concurrent.atomic.AtomicReference;
import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.ApplicationArguments;
import org.springframework.boot.ApplicationRunner;
import org.springframework.hateoas.Resource;
import org.springframework.stereotype.Component;
import uk.ac.ebi.biosamples.client.BioSamplesClient;
import uk.ac.ebi.biosamples.ebeye.gen.*;
import uk.ac.ebi.biosamples.model.Sample;

/*@Component
public class EbEyeBioSamplesDataDumpRunner implements ApplicationRunner {
  private static Logger log = LoggerFactory.getLogger(EbEyeBioSamplesDataDumpRunner.class);
  private static final String BIOSAMPLES = "biosamples";
  private static final String MONGO_SAMPLE = "mongoSample";
  private static final String ENA_LC = "ena";
  private static final String ENA_UC = "ENA";
  private final RefType taxonomyRefType = new RefType();
  @Autowired BioSamplesClient bioSamplesClient;
  @Autowired LoadAttributeSet loadAttributeSet;

  @Value("${spring.data.mongodb.uri}")
  private String mongoUri;

  private Set<String> attributeSet;

  @Override
  public void run(ApplicationArguments args) throws Exception {
    final SimpleDateFormat formatter = new SimpleDateFormat("yyyy-MM-dd", Locale.ENGLISH);
    final DateTimeFormatter dtFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd");
    final MongoClientURI uri = new MongoClientURI(mongoUri);
    final MongoClient mongoClient = new MongoClient(uri);
    final DB db = mongoClient.getDB(BIOSAMPLES);
    final DBCollection coll = db.getCollection(MONGO_SAMPLE);
    final AtomicReference<String> startDate = new AtomicReference<>("");
    final String filePath = "";

    attributeSet = loadAttributeSet.getAllAttributes();

    final List<String> programArguments =
        args.getOptionNames().stream()
            .filter(optionName -> optionName.equals("startDate"))
            .collect(Collectors.toList());

    programArguments.forEach(
        programArgument -> {
          if (programArgument.equals("startDate"))
            startDate.set(args.getOptionValues("startDate").get(0));
        });

    Date startDateFormatted = formatter.parse(String.valueOf(startDate));

    int fileCounter = 1;
    LocalDate startLocalDate = convertToLocalDateViaInstant(startDateFormatted);
    LocalDate endLocalDate = startLocalDate.plusMonths(1).minusDays(1);

    while (endLocalDate.isBefore(LocalDate.now())) {
      File newFile = new File(filePath + "biosd_dump" + fileCounter++ + ".xml");

      log.info(
          "Running for samples with release date starting "
              + startLocalDate.toString()
              + " ending "
              + endLocalDate.toString()
              + " and writing to file "
              + newFile.getPath());

      fetchQueryAndDump(
          coll, startLocalDate.format(dtFormatter), endLocalDate.format(dtFormatter), newFile);

      startLocalDate = startLocalDate.plusMonths(1);
      endLocalDate = startLocalDate.plusMonths(1).minusDays(1);

      // if (fileCounter == 1) break;
    }
  }

  public LocalDate convertToLocalDateViaInstant(Date dateToConvert) {
    return dateToConvert.toInstant().atZone(ZoneId.systemDefault()).toLocalDate();
  }

  private void fetchQueryAndDump(
      final DBCollection coll, final String from, final String until, final File file)
      throws ParseException, JAXBException {
    final List<String> listOfAccessions = getAllDocuments(coll, from, until);

    log.info("Total number of samples to be dumped is : " + listOfAccessions.size());

    List<Sample> samplesList =
        listOfAccessions.stream().map(this::fetchSample).collect(Collectors.toList());

    convertSampleToXml(samplesList, file);
  }

  private static List<String> getAllDocuments(
      final DBCollection col, final String from, final String until) throws ParseException {
    final List<String> listOfAccessions = new ArrayList<>();
    final SimpleDateFormat simpleDateFormat = new SimpleDateFormat("yyyy-MM-dd");
    final Date startDate = simpleDateFormat.parse(from);
    final Date endDate = simpleDateFormat.parse(until);

    final BasicDBObject query =
        new BasicDBObject("release", new BasicDBObject("$gte", startDate).append("$lt", endDate));
    final DBCursor cursor =
        col.find(query).sort(new BasicDBObject("release", OrderBy.ASC.getIntRepresentation()));

    cursor.forEach(
        elem -> {
          final String accession = elem.get("_id").toString();

          log.info("Accession " + accession);
          listOfAccessions.add(elem.get("_id").toString());
        });

    return listOfAccessions;
  }

  public Sample fetchSample(final String accession) {
    Optional<Resource<Sample>> sampleResource = bioSamplesClient.fetchSampleResource(accession);

    return sampleResource.map(Resource::getContent).orElse(null);
  }

  public void convertSampleToXml(final List<Sample> samples, final File file) throws JAXBException {
    DatabaseType databaseType = new DatabaseType();

    databaseType.setName("BioSamples");
    databaseType.setDescription("EBI BioSamples Database");
    databaseType.setEntryCount(samples.size());
    databaseType.setRelease("BioSamples Full Samples Release");
    databaseType.setReleaseDate(new Date().toString());

    AtomicReference<EntriesType> entriesType = new AtomicReference<>(new EntriesType());

    samples.forEach(
        sample -> {
          entriesType.set(getEntries(sample, entriesType.get()));
          databaseType.setEntries(entriesType.get());
        });

    JAXBContext context = JAXBContext.newInstance(DatabaseType.class);

    Marshaller jaxbMarshaller = context.createMarshaller();
    jaxbMarshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);

    jaxbMarshaller.marshal(databaseType, file);
    // jaxbMarshaller.marshal(databaseType, System.out);
  }

  private EntriesType getEntries(final Sample sample, final EntriesType entriesType) {
    final EntryType entryType = new EntryType();

    getEntry(sample, entryType);
    entriesType.getEntry().add(entryType);

    return entriesType;
  }

  private void getEntry(final Sample sample, final EntryType entryType) {
    entryType.setId(sample.getAccession());
    entryType.setName(sample.getName());

    AdditionalFieldsType additionalFieldsType = new AdditionalFieldsType();

    getAdditionalFields(sample, entryType, additionalFieldsType);
    entryType.setAdditionalFields(additionalFieldsType);

    entryType.setDates(getDates(sample));
    entryType.setCrossReferences(getCrossReferences(sample));
  }

  private CrossReferencesType getCrossReferences(final Sample sample) {
    final CrossReferencesType crossReferencesType = new CrossReferencesType();

    sample
        .getExternalReferences()
        .forEach(
            extRef -> {
              RefType refType = new RefType();

              final var url = extRef.getUrl();

              if (url.contains(ENA_LC) || url.contains(ENA_UC)) {
                refType.setDbname(ENA_UC);
                refType.setDbkey(extractEnaAccession(url));
              }

              crossReferencesType.getRef().add(refType);
            });

    crossReferencesType.getRef().add(taxonomyRefType);

    return crossReferencesType;
  }

  private String extractEnaAccession(final String url) {
    return url.substring(36);
  }

  private DatesType getDates(final Sample sample) {
    DatesType datesType = new DatesType();
    DateType dateTypeRelease = new DateType();

    dateTypeRelease.setType("release_date");
    dateTypeRelease.setValue(sample.getReleaseDate());

    DateType dateTypeUpdate = new DateType();

    dateTypeUpdate.setType("update_date");
    dateTypeUpdate.setValue(sample.getUpdateDate());

    datesType.getDate().add(dateTypeRelease);
    datesType.getDate().add(dateTypeUpdate);

    return datesType;
  }

  private void getAdditionalFields(
      final Sample sample,
      final EntryType entryType,
      final AdditionalFieldsType additionalFieldsType) {
    sample
        .getAttributes()
        .forEach(
            attribute -> {
              FieldType fieldType = new FieldType();

              if (attribute.getType().equals("description")) {
                entryType.setDescription(attribute.getValue());
              } else {
                if (sample.getTaxId() != 0) {
                  taxonomyRefType.setDbname("TAXONOMY");
                  taxonomyRefType.setDbkey(sample.getTaxId().toString());
                }

                if (attributeSet.contains(attribute.getType())) {
                  fieldType.setName(
                      removeOtherSpecialCharactersFromAttributeNames(
                          removeSpacesFromAttributeNames(attribute.getType())));
                  fieldType.setValue(attribute.getValue());
                  additionalFieldsType.getFieldOrHierarchicalField().add(fieldType);
                }
              }
            });
  }

  private String removeSpacesFromAttributeNames(String type) {
    return type.trim().replaceAll(" ", "_");
  }

  private String removeOtherSpecialCharactersFromAttributeNames(String type) {
    return type.trim().replaceAll("[^a-zA-Z0-9\\s+_-]", "");
  }
}*/

// One time run for COVID-19 only

@Component
public class EbEyeBioSamplesDataDumpRunner implements ApplicationRunner {
  public static final String ENA_LC = "ena";
  public static final String ENA_UC = "ENA";
  @Autowired BioSamplesClient bioSamplesClient;
  @Autowired EbeyeBioSamplesDataDumpGeneratorDao ebeyeBioSamplesDataDumpGeneratorDao;

  public List<Sample> getSamplesList() {
    Iterable<Resource<Sample>> sampleResources =
        bioSamplesClient.fetchSampleResourceAll("NCBITaxon_2697049");
    List<Sample> sampleList = new ArrayList<>();

    sampleResources.forEach(
        sampleResource -> {
          final Sample sample = sampleResource.getContent();
          /*List<Integer> statusList =
              ebeyeBioSamplesDataDumpGeneratorDao.doGetSampleStatus(sample.getAccession());

          if (statusList != null && statusList.size() > 0) {
            int sampleStatus = statusList.get(0);

            System.out.println(sample.getAccession() + " status " + sampleStatus);

            if (sampleStatus == 5 || sampleStatus == 6) {
              System.out.println(
                  "Sample not added " + sample.getAccession() + " status " + sampleStatus);
            } else {*/
          sampleList.add(sample);
          System.out.println("Sample added " + sample.getAccession());
          /*}
          }*/
        });

    return sampleList;
  }

  public void convertSampleToXml(final List<Sample> samples, final File f) throws JAXBException {
    DatabaseType databaseType = new DatabaseType();

    databaseType.setName("BioSamples");
    databaseType.setDescription("EBI BioSamples Database");
    databaseType.setEntryCount(samples.size());
    databaseType.setRelease("BioSamples COVID-19 Samples Release");
    databaseType.setReleaseDate(new Date().toString());

    AtomicReference<EntriesType> entriesType = new AtomicReference<>(new EntriesType());

    samples.forEach(
        sample -> {
          entriesType.set(getEntries(sample, entriesType.get()));
          databaseType.setEntries(entriesType.get());
        });

    JAXBContext context = JAXBContext.newInstance(DatabaseType.class);

    Marshaller jaxbMarshaller = context.createMarshaller();
    jaxbMarshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);

    jaxbMarshaller.marshal(databaseType, f);
    jaxbMarshaller.marshal(databaseType, System.out);
  }

  private EntriesType getEntries(Sample sample, EntriesType entriesType) {
    EntryType entryType = new EntryType();

    getEntry(sample, entryType);
    entriesType.getEntry().add(entryType);

    return entriesType;
  }

  private void getEntry(Sample sample, EntryType entryType) {
    entryType.setId(sample.getAccession());
    entryType.setName(sample.getName());

    AdditionalFieldsType additionalFieldsType = new AdditionalFieldsType();

    getAdditionalFields(sample, entryType, additionalFieldsType);
    entryType.setAdditionalFields(additionalFieldsType);

    entryType.setDates(getDates(sample));
    entryType.setCrossReferences(getCrossReferences(sample));
  }

  private CrossReferencesType getCrossReferences(Sample sample) {
    CrossReferencesType crossReferencesType = new CrossReferencesType();

    sample
        .getExternalReferences()
        .forEach(
            extRef -> {
              RefType refType = new RefType();

              final String url = extRef.getUrl();

              if (url.contains(ENA_LC) || url.contains(ENA_UC)) {
                refType.setDbname(ENA_UC);
                refType.setDbkey(extractEnaAccession(url));
              }

              crossReferencesType.getRef().add(refType);
            });

    crossReferencesType.getRef().add(getTaxonomyCrossReference(sample.getTaxId()));

    return crossReferencesType;
  }

  private RefType getTaxonomyCrossReference(int taxId) {
    RefType refType = new RefType();

    refType.setDbname("TAXONOMY");
    refType.setDbkey(String.valueOf(taxId));

    return refType;
  }

  private String extractEnaAccession(String url) {
    return url.substring(36);
  }

  private DatesType getDates(Sample sample) {
    DatesType datesType = new DatesType();
    DateType dateTypeRelease = new DateType();

    dateTypeRelease.setType("release_date");
    dateTypeRelease.setValue(sample.getReleaseDate());

    DateType dateTypeUpdate = new DateType();

    dateTypeUpdate.setType("update_date");
    dateTypeUpdate.setValue(sample.getUpdateDate());

    datesType.getDate().add(dateTypeRelease);
    datesType.getDate().add(dateTypeUpdate);

    return datesType;
  }

  private AdditionalFieldsType getAdditionalFields(
      Sample sample, EntryType entryType, AdditionalFieldsType additionalFieldsType) {
    sample
        .getAttributes()
        .forEach(
            attribute -> {
              FieldType fieldType = new FieldType();

              if (attribute.getType().equals("description")) {
                entryType.setDescription(attribute.getValue());
              } else {
                if (attribute.getType().equalsIgnoreCase("host")) {
                  fieldType.setName(
                      removeOtherSpecialCharactersFromAttributeNames(
                          removeSpacesFromAttributeNames("host_scientific_name")));
                  fieldType.setValue(attribute.getValue());
                  additionalFieldsType.getFieldOrHierarchicalField().add(fieldType);
                } else {
                  fieldType.setName(
                      removeOtherSpecialCharactersFromAttributeNames(
                          removeSpacesFromAttributeNames(attribute.getType())));
                  fieldType.setValue(attribute.getValue());
                  additionalFieldsType.getFieldOrHierarchicalField().add(fieldType);
                }
              }
            });

    return additionalFieldsType;
  }

  private String removeSpacesFromAttributeNames(String type) {
    return type.trim().replaceAll(" ", "_");
  }

  private String removeOtherSpecialCharactersFromAttributeNames(String type) {
    return type.trim().replaceAll("[^a-zA-Z0-9\\s+_-]", "");
  }

  @Override
  public void run(ApplicationArguments args) throws Exception {
    File f = new File("biosd-dump_cv_19.xml");
    List<Sample> samplesList = getSamplesList();

    convertSampleToXml(samplesList, f);
  }
}
