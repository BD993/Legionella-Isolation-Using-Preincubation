# Enhancing Culturability of Legionella spp. in Wastewater Samples
# Reproducible analysis. Run from repository root: Rscript Pre-incubation_analysis_v1.0.R
# Requires: ggplot2, dplyr, tidyr, DescTools
# Input: data/*.csv 
# Output: output/* 
# Figure 3 uses Legionella_27_37.csv for both temperature conditions.
# Figure 5 (identification workflow) is intentionally excluded.

required <- c('ggplot2','dplyr','tidyr','DescTools')
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop('Install packages: install.packages(c(',paste(shQuote(missing),collapse=', '),'))')
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr)})
if(!dir.exists('data')) stop('Run from repository root (data/ not found).')
dir.create('output',showWarnings=FALSE,recursive=TRUE)
# Clear previous generated outputs so output/ contains only requested deliverables.
old_outputs <- list.files('output', full.names=TRUE, recursive=FALSE)
old_outputs <- old_outputs[file.info(old_outputs)$isdir == FALSE]
if(length(old_outputs)) unlink(old_outputs)

report <- function(path, title, expr) {
 result <- tryCatch(capture.output(print(force(expr))),
                    error=function(e) paste('ERROR:',conditionMessage(e)))
 cat('\n### ',title,' ###\n',paste(result,collapse='\n'),'\n',
     file=path,append=TRUE,sep='')
 invisible(result)
}
preinc_file <- file.path('output','McNemar_preincubation.txt')
filter_file <- file.path('output','McNemar_filtration.txt')
cochran_file <- file.path('output','Cochran_Q_pretreatments.txt')
writeLines('McNemar tests: paired culture-positive outcomes; chi-squared approximation without continuity correction (correct=FALSE).',preinc_file)
writeLines('McNemar tests: activated sludge, filtered vs. unfiltered; sample-level any-positive across four pre-treatments; correct=FALSE.',filter_file)
writeLines('Cochran Q tests: four unfiltered pre-treatments U, A, T, T+A; additionally all eight treatments for AS.',cochran_file)

read_data <- function(name) {
 path <- file.path('data',name)
 if(!file.exists(path)) stop('Missing input: ',path)
 read.csv2(path,check.names=FALSE,stringsAsFactors=FALSE,na.strings=c('NA','','N/A'))
}
number <- function(x) suppressWarnings(as.numeric(gsub(',','.',as.character(x),fixed=TRUE)))
# ND is a true culture-negative; NA is missing, not a negative.
positive <- function(x) {
 s <- trimws(as.character(x)); out <- rep(NA_integer_,length(s))
 out[!is.na(s) & s=='ND'] <- 0L
 out[!is.na(s) & grepl('^[0-9]+([,.][0-9]+)?$',s)] <- 1L
 out
}
save_jpg <- function(plot,filename,w=8,h=5) {
 ggplot2::ggsave(file.path('output',filename),plot=plot,width=w,height=h,
                units='in',dpi=350,device='jpeg',bg='white',quality=100)
}
theme_pub <- function() theme_classic(base_size=12)+theme(axis.text=element_text(colour='black'),
  strip.background=element_rect(fill='white',colour='black'))

# --- Main Figure 2: culturability, sample-matched tests ---
con <- read_data('Concentrations_Pre-Inc_time.csv')
stopifnot(nrow(con)==50,!anyDuplicated(con$Sample))
con$Type <- ifelse(grepl('_AS$',con$Sample),'AS','TWW')
for(t in c('F','37','LT')) con[[paste0(t,'_positive')]] <- positive(con[[paste0(t,'_Cult')]])
conditions <- c('F','37','LT')
condition_labels <- c(F='Fresh',`37`='Pre-incubated at 37 °C',LT='Pre-incubation at lower temperatures')
rate <- do.call(rbind,lapply(c('All','AS','TWW'),function(group) {
 d <- if(group=='All') con else con[con$Type==group,]
 do.call(rbind,lapply(conditions,function(t) {
  x <- d[[paste0(t,'_positive')]]
  data.frame(Group=group,Condition=condition_labels[[t]],Positive=sum(x==1,na.rm=TRUE),
             N=sum(!is.na(x)),Rate=100*mean(x,na.rm=TRUE))
 }))
}))
rate_all <- subset(rate,Group=='All')
rate_all$Condition <- factor(rate_all$Condition,levels=unname(condition_labels))
# Figure 2 reproduces the manuscript: one blue color, n only, significance bracket.
p2 <- ggplot(rate_all,aes(Condition,Rate))+
 geom_col(width=.30,fill='#0E2B43')+
 geom_text(aes(label=paste0('n = ',N)),vjust=-.75,size=3.8)+
 annotate('segment',x=1,xend=1,y=89,yend=95,linewidth=.35)+
 annotate('segment',x=1,xend=2,y=95,yend=95,linewidth=.35)+
 annotate('segment',x=2,xend=2,y=95,yend=89,linewidth=.35)+
 annotate('text',x=1.5,y=98,label='***',size=6)+
 scale_y_continuous(limits=c(0,105),breaks=seq(0,100,10),
                    labels=function(x)paste0(x,'%'),expand=c(0,0))+
 labs(x=NULL,y='Culturability rate')+theme_classic(base_size=12)+
 theme(axis.text.x=element_text(colour='black',size=10),
       axis.text.y=element_text(colour='black'))
save_jpg(p2,'Figure_2.jpg',9,5)

# McNemar's test: paired complete observations; no continuity correction,
# as in supplied legacy script. Each comparison is separately documented.
mcnemar_pair <- function(d,a,b) {
 x<-d[[paste0(a,'_positive')]]; y<-d[[paste0(b,'_positive')]]
 ok<-complete.cases(x,y)
 tab<-table(factor(x[ok],levels=0:1),factor(y[ok],levels=0:1))
 if(sum(tab[1,2],tab[2,1])==0) return(list(n=sum(ok),table=tab,p.value=NA_real_,note='No discordant pairs'))
 z<-mcnemar.test(tab,correct=FALSE)
 list(n=sum(ok),table=tab,statistic=unname(z$statistic),p.value=z$p.value)
}
for(g in c('All','AS','TWW')) {
 d<-if(g=='All')con else subset(con,Type==g)
 for(pair in list(c('F','37'),c('F','LT')))
  report(preinc_file,paste('McNemar',g,paste(pair,collapse=' vs ')),mcnemar_pair(d,pair[1],pair[2]))
}

# --- Culture pre-treatment comparisons: Cochran's Q ---
tr <- read_data('Cult_pre-treatment.csv')
stopifnot(!anyDuplicated(tr$Sample))
tr$Type <- ifelse(grepl('_AS$',tr$Sample),'AS','TWW')
pretreatments <- c('U','A','T','T.A')
filtered <- c('Uf','Af','Tf','T.Af')
for(stage in c('F','37','LT')) {
 for(group in c('All','AS','TWW')) {
  d <- if(group=='All')tr else tr[tr$Type==group,]
  cols <- paste0(stage,'_',pretreatments)
  mat <- sapply(d[,cols,drop=FALSE],positive)
  mat <- mat[complete.cases(mat),,drop=FALSE]
  report(cochran_file,paste('Cochran Q',stage,group,'unfiltered'),{
   if(nrow(mat)<2) 'Insufficient complete cases' else DescTools::CochranQTest(mat)
  })

 }
 # AS: all eight treatments analysed together (complete cases)
 d <- tr[tr$Type=='AS',]
 mat <- sapply(d[,paste0(stage,'_',c(pretreatments,filtered)),drop=FALSE],positive)
 mat <- mat[complete.cases(mat),,drop=FALSE]
 report(cochran_file,paste('Cochran Q',stage,'AS eight treatments'),
        if(nrow(mat)<2) 'Insufficient complete cases' else DescTools::CochranQTest(mat))
 # Filtration comparison: paired sample-level any-positive, only when all
 # corresponding treatment results are available; NA never becomes zero.
 d <- tr[tr$Type=='AS',]
 u <- sapply(d[,paste0(stage,'_',pretreatments),drop=FALSE],positive)
 f <- sapply(d[,paste0(stage,'_',filtered),drop=FALSE],positive)
 ok <- complete.cases(u)&complete.cases(f)
 if(sum(ok)>0) {
  tab<-table(factor(as.integer(rowSums(u[ok,,drop=FALSE])>0),levels=0:1),
             factor(as.integer(rowSums(f[ok,,drop=FALSE])>0),levels=0:1))
  report(filter_file,paste('AS filtered vs unfiltered',stage),
         if(sum(tab[1,2],tab[2,1])==0) list(n=sum(ok),table=tab,note='No discordance')
         else list(n=sum(ok),table=tab,test=mcnemar.test(tab,correct=FALSE)))
 }
}

# --- Supplementary Figures S2.1–S2.5: temperature profiles ---
# CSV has 240 measurement rows but no timestamps. Therefore the horizontal
# axis is measurement index, not an invented duration in hours or days.
temp <- read_data('T_profiles_LT_37.csv')
# Shared axis ranges across all five supplementary temperature plots.
all_temp_columns <- grep('^[A-E]#[1-5]_(AS|TWW)_27$',names(temp),value=TRUE)
all_temp_values <- number(unlist(temp[,all_temp_columns,drop=FALSE],use.names=FALSE))
common_temp_limits <- c(0, ceiling(max(all_temp_values,na.rm=TRUE)))
common_index_limits <- c(1,nrow(temp))
for(wwtp in LETTERS[1:5]) {
 cols <- grep(paste0('^',wwtp,'#[1-5]_(AS|TWW)_27$'),names(temp),value=TRUE)
 if(!length(cols)) next
 long <- temp[,cols,drop=FALSE] %>% mutate(Measurement=seq_len(n())) %>%
  pivot_longer(-Measurement,names_to='Sample',values_to='Temperature') %>%
  mutate(Temperature=number(Temperature),Event=sub('.*#([1-5])_.*','\\1',Sample)) %>%
  filter(!is.na(Temperature))
 # Profiles for AS and TWW from the same event should coincide; draw one
 # trace per event to avoid double-counting identical incubator readings.
 long <- long %>% group_by(Measurement,Event) %>%
  summarise(Temperature=mean(Temperature),.groups='drop')
 p <- ggplot(long,aes(Measurement,Temperature,colour=Event,group=Event))+
  geom_line(linewidth=.65)+
  scale_colour_manual(values=c('1'='#0072B2','2'='#D55E00','3'='#009E73',
                               '4'='#CC79A7','5'='#E69F00'))+
  coord_cartesian(xlim=common_index_limits,ylim=common_temp_limits)+
  labs(x='Temperature measurement index',y='Temperature (°C)',colour='Sampling',
       title=paste('WWTP',wwtp))+theme_pub()+
   theme(plot.title=element_text(hjust=0.5,face='bold'))
 save_jpg(p,sprintf('Figure_S2.%d.jpg',match(wwtp,LETTERS)),9,5)
}

# --- Main Figure 3: species-by-WWTP temperature plot ---
iso <- read_data('Legionella_27_37.csv')
required_iso <- c('Sample','Average_temperature','Temperature_min','Temperature_max')
stopifnot(all(required_iso %in% names(iso)),!anyDuplicated(iso$Sample))
spcols <- setdiff(names(iso),c(required_iso,'Standard_deviation','No_Legionella'))
iso_long <- iso %>%
 pivot_longer(all_of(spcols),names_to='Species',values_to='Detected') %>%
 filter(!is.na(Detected),tolower(trimws(as.character(Detected)))=='x') %>%
 mutate(WWTP=substr(Sample,1,1),
        Condition=ifelse(grepl('_37$',Sample),'37 °C','Lower temperatures'),
        Average_temperature=number(Average_temperature),
        Temperature_min=number(Temperature_min),
        Temperature_max=number(Temperature_max)) %>%
 filter(is.finite(Average_temperature))

# Explicit species order within each plant; species absent from a plant are omitted.
species_order <- c('L.pneumophila','L.dumoffii','L.rubrilucens/L.erythra',
 'L.londiniensis','Legionella_spp.','uncultured_bacteria','L.feeleii','L.norrlandica')
species_names <- c('L.pneumophila'='L. pneumophila','L.dumoffii'='L. dumoffii',
 'L.rubrilucens/L.erythra'='L. rubrilucens / L. erythra',
 'L.londiniensis'='L. londiniensis','Legionella_spp.'='Unassigned Legionella spp.',
 'uncultured_bacteria'='Uncultured Bacteria','L.feeleii'='L. feeleii',
 'L.norrlandica'='L. norrlandica')
# Plotmath expression labels italicize Latin binomials while keeping generic categories upright.
species_expr <- c('L.pneumophila'='italic("L. pneumophila")',
 'L.dumoffii'='italic("L. dumoffii")',
 'L.rubrilucens/L.erythra'='italic("L. rubrilucens") ~ "/" ~ italic("L. erythra")',
 'L.londiniensis'='italic("L. londiniensis")',
 'Legionella_spp.'='"Unassigned " * italic("Legionella") * " spp."',
 'uncultured_bacteria'='"Uncultured Bacteria"',
 'L.feeleii'='italic("L. feeleii")',
 'L.norrlandica'='italic("L. norrlandica")')
plant_species <- iso_long %>% distinct(WWTP,Species) %>%
 mutate(order=match(Species,species_order)) %>% arrange(WWTP,order) %>%
 group_by(WWTP) %>% mutate(row=row_number()) %>% ungroup()
iso_long <- left_join(iso_long,plant_species,by=c('WWTP','Species'))
# Separate n for lower-temperature and 37 °C observations, matching the source figure.
iso_counts <- iso_long %>% group_by(WWTP,Species,Condition,row) %>%
 summarise(n=n_distinct(Sample),x=ifelse(first(Condition)=='37 °C',37,
            mean(range(Average_temperature))),.groups='drop')
iso_ranges <- iso_long %>% filter(Condition=='Lower temperatures') %>%
 group_by(WWTP,Species,row) %>%
 summarise(xmin=min(Temperature_min,na.rm=TRUE),
           xmax=max(Temperature_max,na.rm=TRUE),.groups='drop')
# Place n below each cross or at the midpoint of each temperature range.
iso_counts <- iso_counts %>% left_join(iso_ranges,by=c('WWTP','Species','row')) %>%
 mutate(x=ifelse(Condition=='Lower temperatures' & is.finite(xmin) & is.finite(xmax),
                 (xmin+xmax)/2,x), y=row+0.33)
# Facet strips are placed to the LEFT of the species-axis labels.
# Free vertical scales keep each WWTP's own set of species.
# Strip backgrounds span the entire WWTP panel (between separators).
plant_species$WWTP <- factor(plant_species$WWTP,levels=LETTERS[1:5],
                            labels=paste('WWTP',LETTERS[1:5]))
iso_long$WWTP <- factor(iso_long$WWTP,levels=LETTERS[1:5],
                        labels=paste('WWTP',LETTERS[1:5]))
iso_ranges$WWTP <- factor(iso_ranges$WWTP,levels=LETTERS[1:5],
                          labels=paste('WWTP',LETTERS[1:5]))
iso_counts$WWTP <- factor(iso_counts$WWTP,levels=LETTERS[1:5],
                          labels=paste('WWTP',LETTERS[1:5]))
# Global unique positions preserve species labels; explicit panel limits align strips and dashed lines.
plant_sizes <- table(plant_species$WWTP)
offsets <- c(0,head(cumsum(as.numeric(plant_sizes)+1),-1))
plant_species$y <- plant_species$row+offsets[as.integer(plant_species$WWTP)]
iso_long <- iso_long %>% left_join(plant_species %>% select(WWTP,Species,y),by=c('WWTP','Species'))
iso_ranges <- iso_ranges %>% left_join(plant_species %>% select(WWTP,Species,y),by=c('WWTP','Species'))
iso_counts <- iso_counts %>% select(-y) %>%
 left_join(plant_species %>% select(WWTP,Species,y),by=c('WWTP','Species')) %>%
 mutate(y=y+0.30)
# Put the count below the cross/line; for the lower-temperature condition
# its x-coordinate is the midpoint of the horizontal temperature range.
# Facet boundaries must be exactly halfway above/below the species rows.
# Without explicit limits, text labels at y + 0.30 expand the panel beyond
# the dashed separator, so the strip boxes and dashed lines don't align.
panel_bounds <- plant_species %>% group_by(WWTP) %>%
 summarise(top=min(y)-0.5,bottom=max(y)+0.5,.groups='drop')
p3 <- ggplot() +
 geom_blank(data=panel_bounds,aes(x=0,y=top))+
 geom_blank(data=panel_bounds,aes(x=40,y=bottom))+
 geom_segment(data=iso_ranges,aes(x=xmin,xend=xmax,y=y,yend=y),linewidth=.65)+
 geom_point(data=iso_long,aes(x=Average_temperature,y=y),shape=4,size=2.5,stroke=.9)+
 geom_text(data=iso_counts,aes(x=x,y=y,label=paste0('n = ',n)),size=3.4,colour='grey35')+
 facet_grid(WWTP~.,scales='free_y',space='free_y',switch='y',drop=FALSE)+
 scale_x_continuous(position='top',limits=c(0,40),breaks=seq(0,40,5),expand=c(0,0))+
 scale_y_reverse(
   breaks=plant_species$y,
   labels=function(b) {
     idx <- match(b,plant_species$y)
     txt <- unname(species_expr[plant_species$Species[idx]])
     txt[is.na(txt)] <- '""'
     as.expression(lapply(txt,function(z) parse(text=z)[[1]]))
   },
   expand=expansion(mult=0,add=0))+
 labs(x='Pre-incubation temperature (°C)',y=NULL)+theme_classic(base_size=13)+
 theme(axis.title.x.top=element_text(face='bold',size=14),
       axis.text.y=element_text(colour='black',size=11),
       axis.text.x.top=element_text(colour='black',size=12),
       axis.line.y=element_line(colour='black'),
       strip.placement='outside',
       strip.text.y.left=element_text(angle=90,face='bold',size=12),
       strip.background=element_rect(fill='grey92',colour='black',linewidth=.4),
       panel.spacing.y=grid::unit(0,'lines'),
       panel.grid=element_blank(),
       plot.margin=margin(12,15,12,8))+
 # Dashed separators coincide with the grey strip's bottom boundary.
 geom_hline(data=panel_bounds %>% filter(as.character(WWTP)!='WWTP E'),
            aes(yintercept=bottom),linetype='dashed',linewidth=.4,
            colour='grey65')
save_jpg(p3,'Figure_3.jpg',10,10)

# --- Main Figure 4: six complete pie charts with left-side condition strips ---
# Counts are transcribed from the manuscript; censored qPCR classification is
# not reconstructed from raw concentrations.
fig4_counts <- data.frame(
 Condition=rep(c('37 °C','Lower temperature conditions'),each=9),
 Group=rep(rep(c('Considering all samples','Considering activated sludge samples',
                 'Considering treated wastewater samples'),each=3),2),
 Change=rep(c('Decrease','Increase','Unchanged'),6),
 Count=c(14,16,6, 8,9,1, 6,7,5, 13,18,5, 7,8,3, 6,10,2))
fig4_counts$Percent <- ave(fig4_counts$Count,
 interaction(fig4_counts$Condition,fig4_counts$Group),FUN=function(z)100*z/sum(z))

fig4_counts$Condition <- factor(fig4_counts$Condition,
 levels=c('37 °C','Lower temperature conditions'))
fig4_counts$Group <- factor(fig4_counts$Group,
 levels=c('Considering all samples','Considering activated sludge samples',
          'Considering treated wastewater samples'))
fig4_counts$Change <- factor(fig4_counts$Change,
 levels=c('Decrease','Increase','Unchanged'))
# Draw the pie wedges as ordinary Cartesian polygons. This avoids clipping
# by coord_polar() inside faceted panels, and creates FULL (not donut) pies.
fig4_counts$RowLabel <- factor(as.character(fig4_counts$Condition),
 levels=levels(fig4_counts$Condition))
fig4_counts$ColLabel <- factor(as.character(fig4_counts$Group),
 levels=levels(fig4_counts$Group))
col_ns <- fig4_counts %>% filter(Condition=='37 °C') %>%
 group_by(Group) %>% summarise(n=sum(Count),.groups='drop')
col_titles <- setNames(paste0(as.character(col_ns$Group),'\nn = ',col_ns$n),
                      as.character(col_ns$Group))
# Calculate angular boundaries separately for each pie.
fig4_wedges <- fig4_counts %>%
 group_by(RowLabel,ColLabel) %>% arrange(Change,.by_group=TRUE) %>%
 mutate(end=2*pi*cumsum(Count)/sum(Count),
        start=lag(end,default=0),mid=(start+end)/2) %>% ungroup()
fig4_polygons <- fig4_wedges %>%
 group_by(RowLabel,ColLabel,Change) %>%
 group_modify(~{
   angles <- seq(.x$start[1],.x$end[1],length.out=120)
   data.frame(x=c(0,cos(angles),0),y=c(0,sin(angles),0))
 }) %>% ungroup()
fig4_labels <- fig4_wedges %>%
 mutate(x=0.67*cos(mid),y=0.67*sin(mid),
        Label=paste0(round(Percent),'%'))
p4 <- ggplot()+
 geom_polygon(data=fig4_polygons,
              aes(x=x,y=y,group=interaction(RowLabel,ColLabel,Change),fill=Change),
              colour='white',linewidth=.4)+
 geom_text(data=fig4_labels,aes(x=x,y=y,label=Label),
           colour='white',size=5)+
 facet_grid(RowLabel~ColLabel,switch='y',labeller=labeller(ColLabel=col_titles))+
 coord_fixed(xlim=c(-1.12,1.12),ylim=c(-1.12,1.12),expand=FALSE)+
 scale_fill_manual(values=c('Decrease'='#126487','Increase'='#C75A0B',
                            'Unchanged'='#00836B'))+
 theme_void(base_size=13)+
 theme(strip.text.x=element_text(face='bold',size=13,margin=margin(b=8)),
       strip.text.y.left=element_text(face='bold',angle=90,size=13),
       strip.placement='outside',strip.background=element_blank(),
       legend.position='bottom', legend.text=element_text(size=13),panel.spacing=grid::unit(0.65,'lines'),
       plot.margin=margin(12,12,12,12))+
 labs(fill=NULL)
save_jpg(p4,'Figure_4.jpg',12,7.8)

cat('Analysis finished. JPG figures and three statistical TXT files are in output/.\n')
