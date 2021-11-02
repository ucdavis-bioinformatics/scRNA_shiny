FROM r-base

# system libraries of general use
## install debian packages
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
libxml2-dev \
libcairo2-dev \
libsqlite3-dev \
libmariadbd-dev \
libpq-dev \
libssh2-1-dev \
unixodbc-dev \
libcurl4-openssl-dev \
libssl-dev

## update system libraries
RUN apt-get update && \
apt-get upgrade -y && \
apt-get clean


## renv.lock file
# COPY /example-app/renv.lock ./renv.lock

# install renv & restore packages
RUN Rscript -e 'install.packages("renv")'
RUN Rscript -e 'install.packages("shiny")'
RUN Rscript -e 'install.packages("ggplot2")'
RUN Rscript -e 'install.packages("Seurat")'
RUN Rscript -e 'install.packages("ggplot2")'
RUN Rscript -e 'install.packages("dplyr")'
RUN Rscript -e 'install.packages("markdown")'
RUN Rscript -e 'install.packages("tidyr")'
RUN Rscript -e 'renv::consent(provided = TRUE)'
RUN Rscript -e 'renv::restore()'
RUN Rscript -e 'memory.limit(size=20000)'


# copy necessary files
## app folder
COPY app.R /app.R
COPY Keller_more_clusterings.rds /Keller_more_clusterings.rds
COPY README.md /README.md


# run app on container start
CMD ["Rscript", "-e", "shiny::runApp('/app.R', host = '0.0.0.0', port = 3838)"]

# expose port
EXPOSE 3838