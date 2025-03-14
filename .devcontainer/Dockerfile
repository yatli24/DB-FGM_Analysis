FROM r-base:4.4.1

# Install system dependencies (adjust as needed)
RUN apt-get update && apt-get install -y \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
  && rm -rf /var/lib/apt/lists/*

# Install required R packages
RUN R -e "install.packages(c('BDgraph','pheatmap','fda','coda', 'igraph', 'jsonlite', 'ggplot2', 'maps', 'mvtnorm', 'fgm', 'ncdf4', 'anytime', 'wavelets', 'doParallel', 'foreach', 'matrixcalc'), repos='http://cran.us.r-project.org')\""

# Set work directory
WORKDIR /home/exouser/DBFGM-Analysis

# Copy your project files into the container
COPY . /home/exouser/DBFGM-Analysis

# Run the simulation study
CMD ["Rscript", "Simulation_study.R"]
