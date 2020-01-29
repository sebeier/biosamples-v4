package uk.ac.ebi.biosamples.ena.amr.service;

import com.fasterxml.jackson.databind.ObjectReader;
import com.fasterxml.jackson.dataformat.csv.CsvMapper;
import com.fasterxml.jackson.dataformat.csv.CsvSchema;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Service;
import uk.ac.ebi.biosamples.client.BioSamplesClient;
import uk.ac.ebi.biosamples.ena.amr.AmrRunner;
import uk.ac.ebi.biosamples.model.Sample;
import uk.ac.ebi.biosamples.model.structured.AbstractData;
import uk.ac.ebi.biosamples.model.structured.amr.AMREntry;
import uk.ac.ebi.biosamples.model.structured.amr.AMRTable;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

@Service
public class EnaAmrDataProcessService {
    private final static Logger log = LoggerFactory.getLogger(EnaAmrDataProcessService.class);

    public void processAmrRows(final List<String> lines, final Sample sample, final BioSamplesClient client) {
        final Set<AbstractData> structuredData = new HashSet<>();
        final AMRTable.Builder amrTableBuilder = new AMRTable.Builder("test");
        /*String[] dilutionMethods = new String[]{"Broth dilution", "Microbroth dilution", "Agar dilution"};
        String[] diffusionMethods = new String[]{"Disc-diffusion", "Neo-sensitabs", "Etest"};*/

        lines.forEach(line -> {
            final CsvMapper mapper = new CsvMapper();
            final CsvSchema schema = mapper.schemaFor(AMREntry.class).withColumnSeparator('\t');
            final ObjectReader r = mapper.readerFor(AMREntry.class).with(schema);

            try {
                AMREntry amrEntry = r.readValue(line);
                amrTableBuilder.addEntry(amrEntry);
            } catch (IOException e) {
                e.printStackTrace();
            }
        });

        structuredData.add(amrTableBuilder.build());
        Sample sampleNew = Sample.Builder.fromSample(sample).withData(structuredData).build();
        client.persistSampleResource(sampleNew);

        log.info("Submitted sample " + sample.getAccession() + " with structured data");
    }

    private String removeBioSampleId(String line) {
        return line.substring(line.indexOf(AmrRunner.TAB));
    }

    public List<String> processAmrLines(BufferedReader bufferedReader) {
        return bufferedReader.lines().skip(1).map(this::removeBioSampleId).collect(Collectors.toList());
    }
}