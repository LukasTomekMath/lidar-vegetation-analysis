library(readr)
library(dplyr)
library(tidyr)

#NACITANIE
biotopy_nazvy <- read_csv("C:/Users/Maria/Desktop/dilinovka/databaza/biotopy_nazvy.csv")
biotopy_segmentov <- read_csv("C:/Users/Maria/Desktop/dilinovka/databaza/biotopy_segmentov.csv")
lokality_id <- read_csv("C:/Users/Maria/Desktop/dilinovka/databaza/lokality_id.csv")


#odtial chcem len relevantne, vela stlpcov tam je inak
lokality_id <- lokality_id %>%
  select(all_of(c("LokalitaID","LokalitaNazov","Hectare","Area","Perimeter","OrografickeCelky")))
 
#chcem info z GisPrvokID_kb separante a zmazat tie co nemaju kod pralesa
biotopy_segmentov <- biotopy_segmentov %>%
  separate_wider_delim(
    cols = GisPrvokID_kb,
    delim = "_",
    names = c("LokalitaID", "SegmentID", "Datum", "MapovatelHlavnyID", "KodPralesovejLokality"),
    too_few = "align_start"
  ) %>%
  filter(!is.na(KodPralesovejLokality)) 

#chcem spojit biotopy segmentov aj s nazvami iotopov nech to vidim
spolu <- biotopy_segmentov %>%
  left_join(biotopy_nazvy, by = c("BiotopID" = "BiotopID")) 

#este to chcem spojit aj s lokalitaID
spolu <- spolu %>%
mutate(LokalitaID = as.character(LokalitaID)) %>% 
  left_join(
    lokality_id %>% mutate(LokalitaID = as.character(LokalitaID)), 
    by = "LokalitaID"
  ) %>%
  mutate(LokalitaID = as.numeric(LokalitaID))

spolu<- spolu %>% 
  select(-c(Datum,MapovatelHlavnyID,Hectare,Area,Perimeter,zapisID_kb))

write_excel_csv(spolu, "C:/Users/Maria/Desktop/dilinovka/databaza/spolu_export.csv")
