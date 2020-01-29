package uk.ac.ebi.biosamples.ena.amr;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.ApplicationArguments;
import org.springframework.boot.ApplicationRunner;
import org.springframework.hateoas.Resource;
import org.springframework.stereotype.Component;
import uk.ac.ebi.biosamples.PipelinesProperties;
import uk.ac.ebi.biosamples.client.BioSamplesClient;
import uk.ac.ebi.biosamples.ena.amr.service.EnaAmrDataProcessService;
import uk.ac.ebi.biosamples.model.Sample;
import uk.ac.ebi.biosamples.utils.AdaptiveThreadPoolExecutor;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.Future;

@Component
public class AmrRunner implements ApplicationRunner {
    private static final String BSD_SAMPLE_PREFIX = "SA";
    private static final String FTP = "ftp";
    private static final String MD_5 = "md5";
    private static final String HTTP = "http://";
    public static final String TAB = "\t";
    private static final String antibiogram = "\"AMR_ANTIBIOGRAM\"";
    private static final String URL = "https://www.ebi.ac.uk/ena/portal/api/search?result=analysis&query=analysis_type=" + antibiogram + "&dataPortal=pathogen&dccDataOnly=false&fields=analysis_accession,country,region,scientific_name,location,sample_accession,tax_id,submitted_ftp,first_public,last_updated&sortFields=scientific_name,country&limit=0";
    private final static Logger log = LoggerFactory.getLogger(AmrRunner.class);
    private final static Map<String, Future<Void>> futures = new LinkedHashMap<>();

    @Autowired
    EnaAmrDataProcessService enaAmrDataProcessService;
    @Autowired
    BioSamplesClient bioSamplesClient;


    @Override
    public void run(final ApplicationArguments args) {
        List<AccessionFtpUrlPair> pairList = new ArrayList<>();

        try {
            pairList = requestHttpAndGetAccessionFtpUrlPairs();
        } catch (Exception e) {
            log.info("An exception occured while fetching AMR data from ENA API " + e.getMessage());
        }

        downloadFtpContent(pairList);
    }

    private static List<AccessionFtpUrlPair> requestHttpAndGetAccessionFtpUrlPairs() throws Exception {
        final URL enaApiUrl = new URL(AmrRunner.URL);
        final HttpURLConnection conn = (HttpURLConnection) enaApiUrl.openConnection();
        List<AccessionFtpUrlPair> pairList = new ArrayList<>();

        try {
            if (getResponseFromEnaApi(conn) == 200) {
                pairList = doGetAccessionFtpUrlPairs(enaApiUrl);
            }
        } catch (final Exception e) {
            throw new RuntimeException(e);
        } finally {
            conn.disconnect();
        }

        return pairList;
    }

    private static int getResponseFromEnaApi(final HttpURLConnection conn) throws IOException {
        int response;

        conn.setRequestMethod("GET");
        conn.connect();
        response = conn.getResponseCode();

        return response;
    }

    private static List<AccessionFtpUrlPair> doGetAccessionFtpUrlPairs(final URL url) {
        final List<AccessionFtpUrlPair> accessionFtpUrlPairs = new ArrayList<>();

        try {
            BufferedReader bufferedReader = getReader(url);

            bufferedReader.lines().forEach(line -> accessionFtpUrlPairs.add(getAccessionFtpUrlPair(line)));
        } catch (final IOException e) {
            log.info("Failed to get and parse accession and FTP pairs for URL " + url.toString());
        }

        return accessionFtpUrlPairs;
    }

    private static BufferedReader getReader(final URL url) throws IOException {
        return new BufferedReader(new InputStreamReader(url.openConnection().getInputStream()));
    }

    private static AccessionFtpUrlPair getAccessionFtpUrlPair(final String line) {
        final StringTokenizer tokenizer = new StringTokenizer(line, TAB);
        final AccessionFtpUrlPair accessionFtpUrlPair = new AccessionFtpUrlPair();

        while (tokenizer.hasMoreTokens()) {
            final String value = tokenizer.nextToken();

            if (value.startsWith(BSD_SAMPLE_PREFIX)) {
                accessionFtpUrlPair.setAccession(value);
            }

            if (value.startsWith(FTP)) {
                dealWithSemicolon(value, accessionFtpUrlPair);
            }
        }

        return accessionFtpUrlPair;
    }

    private static void dealWithSemicolon(final String value, final AccessionFtpUrlPair accessionFtpUrlPair) {
        final int index = value.indexOf(';');
        final String option1 = value.substring(index + 1);
        final String option2 = value.substring(0, index);

        if (!option1.endsWith(MD_5)) {
            accessionFtpUrlPair.setFtpUrl(HTTP + option1);
        } else {
            accessionFtpUrlPair.setFtpUrl(HTTP + option2);
        }
    }

    private void downloadFtpContent(final List<AccessionFtpUrlPair> pairList) {
        try (final AdaptiveThreadPoolExecutor executorService = AdaptiveThreadPoolExecutor.create(100, 10000, true,
                1, 10)) {

            pairList.forEach(pair -> {
                try {
                    executorService.submit(new EnaAmrCallable(new URL(pair.getFtpUrl()), pair.getAccession()));
                } catch (MalformedURLException e) {
                    log.info("FTP URL not correctly formed " + pair.getFtpUrl());
                }
            });
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    class EnaAmrCallable implements Callable<Void> {
        URL url;
        String accession;

        EnaAmrCallable(final URL url, final String accession) {
            this.url = url;
            this.accession = accession;
        }

        @Override
        public Void call() {
            return fetchSampleAndProcessAmrData(url, accession);
        }

        private Void fetchSampleAndProcessAmrData(final URL url, final String accession) {
            try {
                final Optional<Resource<Sample>> sample = bioSamplesClient.fetchSampleResource(accession, Optional.of(new ArrayList<>()));

                if (sample.isPresent()) {
                    enaAmrDataProcessService.processAmrRows(enaAmrDataProcessService.processAmrLines(getReader(url)), sample.get().getContent(), bioSamplesClient);
                } else {
                    log.info(accession + " doesn't exist");
                }
            } catch (final IOException ioe) {
                log.info("Couldn't process AMR data for " + accession);
            }

            return null;
        }
    }
}