# =======================================================================
# Builds data/births.rda from the INEC live birth file ENV_2024.csv
# (Estadisticas de Nacidos Vivos 2024, INEC Ecuador).
# Run from the package root. ENV_2024.csv is not part of the package.
# =======================================================================

library(readr)
library(dplyr)

env <- read_delim("ENV_2024.csv", delim = ";", locale = locale(encoding = "UTF-8"))

dat <- env %>%
  mutate(
    fecha_nac  = as.Date(fecha_nac,  format = "%m/%d/%Y"),
    fecha_insc = as.Date(fecha_insc, format = "%Y/%m/%d"),
    lag        = as.numeric(fecha_insc - fecha_nac),
    home       = ifelse(lugar_ocur == "Casa", 1, 0),
    rural      = ifelse(area_res == "Rural", 1, 0),
    indigena   = ifelse(etnia == "Ind\u00edgena", 1, 0),
    out_prov   = ifelse(prov_insc != prov_nac, 1, 0)
  ) %>%
  filter(
    !is.na(lag),
    lag >= 0,
    !lugar_ocur %in% c("Sin informacion", "Sin informaci\u00f3n")
  )

# Jan 2024 births
coh <- dat %>%
  filter(anio_nac == 2024, as.integer(format(fecha_nac, "%m")) <= 1) %>%
  filter(!is.na(home), !is.na(rural), !is.na(indigena), !is.na(out_prov))

# keep only the variables used in the paper, with English names
births <- data.frame(
  birth_date        = coh$fecha_nac,
  registration_date = coh$fecha_insc,
  lag               = coh$lag,
  home              = coh$home,
  rural             = coh$rural,
  indigenous        = coh$indigena,
  out_province      = coh$out_prov
)

save(births, file = "data/births.rda", compress = "xz")
