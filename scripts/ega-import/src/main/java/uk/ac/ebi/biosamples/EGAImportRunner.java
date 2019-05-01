package uk.ac.ebi.biosamples;

import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.ApplicationArguments;
import org.springframework.boot.ApplicationRunner;
import org.springframework.hateoas.Resource;
import org.springframework.stereotype.Component;
import uk.ac.ebi.biosamples.client.BioSamplesClient;
import uk.ac.ebi.biosamples.model.Attribute;
import uk.ac.ebi.biosamples.model.ExternalReference;
import uk.ac.ebi.biosamples.model.Sample;
import uk.ac.ebi.biosamples.ols.OlsProcessor;
import uk.ac.ebi.biosamples.ols.OlsResult;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.time.Instant;
import java.util.*;
import java.util.stream.Collectors;

@Component
public class EGAImportRunner implements ApplicationRunner {
    private static final Logger LOG = LoggerFactory.getLogger(EGAImportRunner.class);
    private static final String EGA_DATASET_BASE_URL = "https://ega-archive.org/datasets/";
    private static final String EGA_SAMPLE_BASE_URL = "https://ega-archive.org/metadata/v2/samples/";
    private static final Set<String> UNKNOWN_TERMS = new HashSet<>(Arrays.asList("n/a", "na", "n.a", "none",
            "unknown", "--", ".", "null", "missing", "[not reported]", "[not requested]", "not applicable",
            "not_applicable", "not collected", "not specified", "not known", "not reported", "missing: not provided"));
    private static final String ATTRIBUTE_PHENOTYPE = "phenotype";
    private static final String ATTRIBUTE_SEX = "sex";

    private final Attribute organism;
    private final BioSamplesClient bioSamplesClient;
    private final OlsProcessor olsProcessor;
    private final Map<String, Optional<OlsResult>> olsResultCache;

    @Autowired
    public EGAImportRunner(BioSamplesClient bioSamplesClient, OlsProcessor olsProcessor) {
        this.bioSamplesClient = bioSamplesClient;
        this.olsProcessor = olsProcessor;

        olsResultCache = new HashMap<>();
        organism = Attribute.build("organism", "Homo sapiens", "http://purl.obolibrary.org/obo/NCBITaxon_9606", null);
    }

    @Override
    public void run(ApplicationArguments args) {
        if (args.getSourceArgs().length < 1) {
            LOG.error("Please specify a data folder as a program argument");
            throw new IllegalArgumentException("Please specify a data folder as a program argument");
        }

        final String dataFolderUrl = args.getSourceArgs()[0];
        final String datasetDuoUrl = dataFolderUrl + "1.datasets_duo.csv";
        final String sampleDataUrl = dataFolderUrl + "1.2.sanger_released_samples.csv";
        final String phenotypeIriFile = dataFolderUrl + "sanger_datasets_public_phenotype_hpo.csv";

        Map<String, SortedSet<String>> datasetToDuoCodesMap = loadDuoCodeMap(datasetDuoUrl);
//        Map<String, List<OlsResult>> phenotypeIriMap = loadPhenotypeIriMap(phenotypeIriFile);
        Map<String, List<OlsResult>> phenotypeIriMap = new HashMap<>();//todo remove this and uncomment above

        try {
            String[] headers = {"biosample_id", "sample_id", "dataset_id", "phenotype", "gender"};
            FileReader in = new FileReader(sampleDataUrl);
            Iterable<CSVRecord> records = CSVFormat.DEFAULT
                    .withHeader(headers)
                    .withFirstRecordAsHeader()
                    .parse(in);
            for (CSVRecord record : records) {
                String accession = record.get("biosample_id");
                String egaId = record.get("sample_id");
                String datasetId = record.get("dataset_id");
                String phenotype = record.get("phenotype");
                String sex = record.get("sex");
                //String sex = "unknown";

                SortedSet<String> duoCodes = datasetToDuoCodesMap.get(datasetId);
                List<OlsResult> phenotypeIris = phenotypeIriMap.get(phenotype);

                processSampleRecord(accession, egaId, datasetId, phenotype, sex, duoCodes, phenotypeIris);
            }
        } catch (FileNotFoundException e) {
            LOG.error("Couldn't read file: " + datasetDuoUrl, e);
        } catch (IOException e) {
            LOG.error("Failed to parse the file: " + datasetDuoUrl, e);
        }
    }

    private void processSampleRecord(String accession, String egaId, String datasetId, String phenotype, String sex,
                                     SortedSet<String> duoCodes, List<OlsResult> phenotypeIris)
            throws JsonProcessingException {

        final ObjectMapper jsonMapper = new ObjectMapper();
        Optional<Resource<Sample>> sampleResource = bioSamplesClient.fetchSampleResource(accession);
        if (sampleResource.isPresent()) {
            Sample sample = sampleResource.get().getContent();
            LOG.info("Original sample: {}", jsonMapper.writeValueAsString(sample));
            if (sample.getAttributes().size() != 2) {
                LOG.warn("{}: Attributes size != 2, Attributes {}", sample.getAccession(), sample.getAttributes());
            }

            //remove extra attributes from migration (deleted and other-migrated from....)
            removeMigrationRelatedAttributes(sample);

            Sample.Builder sampleBuilder = Sample.Builder.fromSample(sample)
                    .addAttribute(Attribute.build("ega dataset id", datasetId))
                    .addAttribute(Attribute.build("ega sample id", egaId))
                    .addAttribute(organism)
                    .addExternalReference(ExternalReference.build(EGA_DATASET_BASE_URL + datasetId, duoCodes))
                    .addExternalReference(ExternalReference.build(EGA_SAMPLE_BASE_URL + egaId))
//                    .withRelease(Instant.now());
                    .withRelease(Instant.parse("2019-04-01T00:00:01Z"));

            //ignore unknown, n/a terms
            if (phenotype == null || "".equals(phenotype) || UNKNOWN_TERMS.contains(phenotype.toLowerCase())) {
                LOG.info("Ignoring phenotype as it contains {}", phenotype);
            } else {
                Attribute attributePhenotype = populateAttribute(phenotype, phenotypeIris, ATTRIBUTE_PHENOTYPE);
                sampleBuilder.addAttribute(attributePhenotype);
            }
            if (sex == null || "".equals(sex) || UNKNOWN_TERMS.contains(sex.toLowerCase())) {
                LOG.info("Ignoring sex as it contains {}", sex);
            } else {
                Attribute attributeSex = populateAttribute(sex.toLowerCase(), getSexOntology(sex), ATTRIBUTE_SEX);
                sampleBuilder.addAttribute(attributeSex);
            }

            Sample updatedSample = sampleBuilder.build();
            bioSamplesClient.persistSampleResource(updatedSample);
            LOG.info("Updated sample: {}", jsonMapper.writeValueAsString(updatedSample));
        } else {
            LOG.warn("Sample not found in biosamples: {}", accession);
        }
    }

    private Map<String, SortedSet<String>> loadDuoCodeMap(String datasetDuoUrl) {
        Map<String, SortedSet<String>> datasetToDuoCodesMap = new HashMap<>();

        try {
            String[] headers = {"dataset_id", "duo_codes_json"};
            FileReader in = new FileReader(datasetDuoUrl);
            Iterable<CSVRecord> records = CSVFormat.DEFAULT
                    .withHeader(headers)
                    .withFirstRecordAsHeader()
                    .parse(in);
            for (CSVRecord record : records) {
                String datasetId = record.get("dataset_id");
                String duoCodesArrayAsString = record.get("duo_codes_json");
                String[] duoCodes = duoCodesArrayAsString.replaceAll("[\"\\[\\] ]", "").split(",");

                datasetToDuoCodesMap.put(datasetId,
                        new TreeSet<>(Arrays.stream(duoCodes).map(s -> "DUO:" + s).collect(Collectors.toList())));
            }
        } catch (FileNotFoundException e) {
            LOG.error("couldn't read file: " + datasetToDuoCodesMap, e);
        } catch (IOException e) {
            LOG.error("Failed to parse the file: " + datasetToDuoCodesMap, e);
        }

        return datasetToDuoCodesMap;
    }

    private Map<String, List<OlsResult>> loadPhenotypeIriMap(String phenotypeIriFile) {
        Map<String, List<OlsResult>> phenotypeIriMap = new HashMap<>();

        try {
            String[] headers = {"public_phenotype", "mapped_term", "hpo_id", "efo_id"};
            FileReader in = new FileReader(phenotypeIriFile);
            Iterable<CSVRecord> records = CSVFormat.DEFAULT
                    .withHeader(headers)
                    .withFirstRecordAsHeader()
                    .parse(in);
            for (CSVRecord record : records) {
                String publicPhenotype = record.get("public_phenotype");
                String mappedPhenotype = record.get("mapped_term");
                String hpoId = record.get("hpo_id");
                String efoId = record.get("efo_id");
                List<OlsResult> iriSet = new ArrayList<>();

                if (hpoId != null && !"".equals(hpoId)) {
                    Optional<OlsResult> olsResult = getOlsMappedTerm(hpoId);
                    olsResult.ifPresent(iriSet::add);
                }
                if (efoId != null && !"".equals(efoId)) {
                    Optional<OlsResult> olsResult = getOlsMappedTerm(efoId);
                    olsResult.ifPresent(iriSet::add);
                }

                phenotypeIriMap.put(publicPhenotype, iriSet);
            }
        } catch (FileNotFoundException e) {
            LOG.error("couldn't read file: " + phenotypeIriFile, e);
        } catch (IOException e) {
            LOG.error("Failed to parse the file: " + phenotypeIriFile, e);
        }

        return phenotypeIriMap;
    }

    private void removeMigrationRelatedAttributes(Sample sample) {
        List<Attribute> attributesToRemove = new ArrayList<>();
        for (Attribute attribute : sample.getAttributes()) {
            if (attribute.getType().equals("deleted") ||
                    (attribute.getType().equals("other") && attribute.getValue().startsWith("migrated from"))) {
                attributesToRemove.add(attribute);
                LOG.info("Removing attribute {}={} from original sample: {}", attribute.getType(), attribute.getValue(), sample.getAccession());
            } else if (attribute.getType().equals("phenotype")) {
                attributesToRemove.add(attribute);
                LOG.warn("Removing attribute phenotype={} from original sample: {}", attribute.getValue(), sample.getAccession());
            } else if (attribute.getType().equals("organism")) {
                attributesToRemove.add(attribute);
                LOG.warn("Removing attribute organism={} from original sample: {}", attribute.getValue(), sample.getAccession());
            } else if (attribute.getType().equals("sex")) {
                attributesToRemove.add(attribute);
                LOG.warn("Removing attribute sex={} from original sample: {}", attribute.getValue(), sample.getAccession());
            }//todo ega dataset and sample.
        }
        for (Attribute attribute : attributesToRemove) {
            sample.getAttributes().remove(attribute);
        }
    }

    private Attribute populateAttribute(String phenotype, List<OlsResult> attributeIris, String attributeType) {
        Optional<OlsResult> olsMappedTerm = getOlsMappedTerm(phenotype);
        Attribute attribute;

        List<String> iris = new ArrayList<>();
        if (attributeIris != null && !attributeIris.isEmpty()) {
            for (OlsResult o : attributeIris) {
                iris.add(o.getIri());
            }
        }

        if (olsMappedTerm.isPresent()) {
            iris.add(olsMappedTerm.get().getIri());
            attribute = Attribute.build(attributeType, olsMappedTerm.get().getLabel(), iris, null);
        } else {
            attribute = Attribute.build(attributeType, phenotype, iris, null);
        }

        return attribute;
    }

    private Optional<OlsResult> getOlsMappedTerm(String termToMap) {
        Optional<OlsResult> olsMappedTerm = Optional.empty();
        if (termToMap.matches("^[A-Za-z]+[_:\\-][0-9]+$")) {
            if (olsResultCache.containsKey(termToMap)) {
                olsMappedTerm = olsResultCache.get(termToMap);
            } else {
                olsMappedTerm = olsProcessor.queryForOlsObject(termToMap);
                olsResultCache.put(termToMap, olsMappedTerm);
            }

            if (olsMappedTerm.isPresent()) {
                LOG.info("OLS mapped term {} into {}", termToMap, olsMappedTerm.get().getIri());
            } else {
                LOG.warn("Could not find term({}) in OLS", termToMap);
            }
        }

        return olsMappedTerm;
    }

    private List<OlsResult> getSexOntology(String sex) {
        List<OlsResult> olsResults;
        switch (sex.toLowerCase()) {
            case "male":
                olsResults = Collections.singletonList(new OlsResult("male", "http://purl.obolibrary.org/obo/PATO_0000384"));
                break;
            case "female":
                olsResults = Collections.singletonList(new OlsResult("female", "http://purl.obolibrary.org/obo/PATO_0000383"));
                break;
            default:
                olsResults = null;
                break;
        }
        return olsResults;
    }
}
