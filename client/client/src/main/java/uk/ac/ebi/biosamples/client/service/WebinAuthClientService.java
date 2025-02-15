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
package uk.ac.ebi.biosamples.client.service;

import com.auth0.jwt.JWT;
import com.auth0.jwt.interfaces.DecodedJWT;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.net.URI;
import java.util.Date;
import java.util.List;
import java.util.Optional;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.web.client.RestTemplateBuilder;
import org.springframework.http.*;
import org.springframework.web.client.RestTemplate;
import uk.ac.ebi.biosamples.client.model.auth.AuthRealm;
import uk.ac.ebi.biosamples.client.model.auth.AuthRequestWebin;

public class WebinAuthClientService implements ClientService {
  private Logger log = LoggerFactory.getLogger(getClass());

  private final RestTemplate restTemplate;

  private final URI webinAuthUri;
  private final String username;
  private final String password;
  private final List<AuthRealm> authRealms;

  private Optional<String> jwt = Optional.empty();
  private Optional<Date> expiry = Optional.empty();

  @Autowired ObjectMapper objectMapper;

  public WebinAuthClientService(
      RestTemplateBuilder restTemplateBuilder,
      URI webinAuthUri,
      String username,
      String password,
      List<AuthRealm> authRealms) {
    this.restTemplate = restTemplateBuilder.build();
    this.webinAuthUri = webinAuthUri;
    this.username = username;
    this.password = password;
    this.authRealms = authRealms;
  }

  public synchronized String getJwt() {
    if (username == null
        || username.trim().length() == 0
        || password == null
        || password.trim().length() == 0
        || authRealms.size() == 0) {
      return null;
    }

    if (!jwt.isPresent() || (expiry.isPresent() && expiry.get().before(new Date()))) {
      final AuthRequestWebin authRequestWebin =
          new AuthRequestWebin(username, password, authRealms);
      HttpHeaders headers = new HttpHeaders();
      headers.setContentType(MediaType.APPLICATION_JSON);
      HttpEntity<String> entity;

      try {
        entity = new HttpEntity<>(objectMapper.writeValueAsString(authRequestWebin), headers);

        ResponseEntity<String> responseEntity =
            restTemplate.exchange(webinAuthUri, HttpMethod.POST, entity, String.class);

        final HttpStatus statusCode = responseEntity.getStatusCode();

        if (statusCode == HttpStatus.OK || statusCode == HttpStatus.CREATED) {
          jwt = Optional.of(responseEntity.getBody());
        } else jwt = Optional.empty();

        DecodedJWT decodedJwt = JWT.decode(jwt.orElse(null));
        expiry = Optional.of(decodedJwt.getExpiresAt());
      } catch (Exception e) {
        throw new RuntimeException(e);
      }

      log.info("jwt = " + jwt);
    }

    return jwt.get();
  }

  @Override
  public boolean isWebin() {
    return true;
  }
}
