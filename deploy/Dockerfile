FROM adviselipidomics_base
COPY renv.lock.prod renv.lock
RUN R -e 'renv::restore()'
COPY ADViSELipidomics_*.tar.gz /app.tar.gz
RUN R -e 'remotes::install_local("/app.tar.gz",upgrade="never")'
RUN rm /app.tar.gz
EXPOSE 80
CMD R -e "options('shiny.port'=80,shiny.host='0.0.0.0');library(ADViSELipidomics);ADViSELipidomics::run_ADViSELipidomics()"
