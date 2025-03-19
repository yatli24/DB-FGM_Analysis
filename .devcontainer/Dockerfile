FROM r-base:4.4.1

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
  && rm -rf /var/lib/apt/lists/*

# Install required R packages
RUN R -e "install.packages(c('BDgraph','pheatmap','fda','coda','igraph','gt','jsonlite','ggplot2','maps','mvtnorm','fgm','ncdf4','anytime','wavelets','doParallel','foreach','matrixcalc'), repos='http://cran.us.r-project.org')"

# Set work directory
WORKDIR /home/exouser/DB-FGM_Analysis

# Keep the container running
CMD ["tail", "-f", "/dev/null"]
