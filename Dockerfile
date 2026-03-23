FROM rocker/shiny:latest

# Install dependencies (including libpng/libjpeg for Seurat spatial image rendering)
RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    libglpk-dev \
    libpng-dev \
    libjpeg-dev \
    python3 python3-pip python3-venv \
    vim nano \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Create and activate Python env if needed
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN pip install pandas

# Install R packages — wrapper ensures the build fails loudly if any package
# is missing after install (install.packages() only warns on failure by default)
RUN R -e '
  install_required <- function(pkgs) {
    install.packages(pkgs)
    missing <- pkgs[!pkgs %in% rownames(installed.packages())]
    if (length(missing)) stop("Failed to install: ", paste(missing, collapse = ", "))
    invisible(NULL)
  }
  install_required(c(
    "plotly", "igraph", "Matrix", "DT", "ggplot2", "dplyr", "stringr",
    "shiny", "shinymanager", "shinyWidgets", "bslib", "shinycssloaders",
    "bsicons", "periscope2", "shinyjs", "tidyr", "DBI", "RSQLite",
    "jsonlite", "png"
  ))
'
# Install Seurat ecosystem separately (large; keep in its own layer)
RUN R -e '
  install.packages(c("SeuratObject", "Seurat"), repos="https://cloud.r-project.org/")
  missing <- c("SeuratObject", "Seurat")[!c("SeuratObject", "Seurat") %in% rownames(installed.packages())]
  if (length(missing)) stop("Failed to install: ", paste(missing, collapse = ", "))
'
RUN R -e '
  install.packages("VAM")
  if (!"VAM" %in% rownames(installed.packages())) stop("Failed to install VAM")
'

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
