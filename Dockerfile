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
RUN R -e 'install.packages(c("plotly", "igraph", "Seurat", "Matrix", "SeuratObject", "DT", "ggplot2", "dplyr", "stringr", "shiny", "shinymanager", "shinyWidgets", "bslib", "shinycssloaders", "bsicons", "VAM", "periscope2", "shinyjs", "tidyr"))'

# Copy Shiny app code
COPY ./app_code/ /srv/shiny-server/atlas/
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf
RUN mkdir -p /var/log/shiny-server && \
    chown -R shiny:shiny /var/log/shiny-server && \
    mkdir -p /srv/shiny-server/atlas/app_cache && \
    chmod 777 /srv/shiny-server/atlas/app_cache && \
    mkdir -p /srv/shiny-server/atlas/data && \
    chown -R shiny:shiny /srv/shiny-server/atlas/data && \
    chmod 755 /srv/shiny-server/atlas/data && \
    find /srv/shiny-server/atlas/data -type f -name "*.sqlite*" -exec chmod 664 {} \; && \
    find /srv/shiny-server/atlas/data -type f -name "*.db" -exec chmod 664 {} \;


# Expose Shiny port
EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
