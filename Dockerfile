FROM rocker/shiny:latest

# Install dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    libglpk-dev \
    python3 python3-pip python3-venv \
    vim nano

# Create and activate Python env if needed
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN pip install pandas

# Install R packages
RUN R -e 'install.packages(c("plotly", "igraph", "Seurat", "Matrix", "SeuratObject", "DT", "ggplot2", "dplyr", "stringr", "shiny", "shinyWidgets", "bslib", "shinycssloaders", "bsicons", "VAM", "periscope2", "shinyjs"))'

# Copy Shiny app code
COPY ./app_code/ /srv/shiny-server/atlas/
COPY ./imgs/ /srv/shiny-server/atlas/imgs/
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf
RUN mkdir -p /var/log/shiny-server && \
    chown -R shiny:shiny /var/log/shiny-server
    mkdir -p /srv/shiny-server/atlas/app_cache
    chmod 777 /srv/shiny-server/atlas/app_cache

# Expose Shiny port
EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
