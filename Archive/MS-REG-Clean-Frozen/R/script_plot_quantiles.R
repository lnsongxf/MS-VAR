if (!require("devtools")) install.packages("devtools")
library(devtools)
library(httr)

library(tidyverse)
library(ggstream)
library(patchwork)
library(wesanderson)
library(tinter)
library(ggthemes)
library(latex2exp)

set_config(use_proxy('wwwproxy.frb.gov',8080), override = FALSE)
set_config(user_agent('Lynx'), override = FALSE)


# Palette
pal2 <- wes_palette("FantasticFox1",5)
pal <- wes_palette("Chevalier1")[c(4, 1:3)]

model = "Direct"
# Read in data
df <- readr::read_csv(paste0('dataquantiles',model,'.csv'))
df[["dates"]] <- as.Date(df[["dates"]], "%d-%b-%Y")

# Fonts
f1 = "Piazzolla SC"

t2009<-subset(df,dates==as.Date(c("2007-01-01")))
t2020<-subset(df,dates==as.Date(c("2020-01-01")))

p3 <-ggplot(df) +
  geom_ribbon(aes(x=dates, ymax = dYsim_90, ymin = dYsim_10), alpha = 0.2,
              fill = pal[3], color = "transparent")+
  geom_line(data=df, aes(x=dates, y=dYsim_90_cf1), linetype='dashed',color="steelblue4", lwd = 2.0) +
  geom_line(data=df, aes(x=dates, y=dYsim_10_cf1), linetype='dashed',color=pal2[5], lwd = 2.0) +
  geom_line(aes(dates, dYsim_90), size = 0.5,color=colorspace::lighten(pal[3], 0.5)) +  
  geom_line(aes(dates, dYsim_10), size = 0.5,color=colorspace::lighten(pal[3], 0.5)) +
  annotate("text", x = t2009$dates, y = 0.5*(t2009$dYsim_10+t2009$dYsim_90), label = "Estimated", size = 3,color=colorspace::darken(pal[3],0.1),hjust=0)+
  annotate("text", x = t2009$dates, y = t2009$dYsim_10_cf1*0.70, label = "Q10: no MF", size = 3,color=colorspace::darken(pal2[5],0.1),hjust=0,vjust=1)+
  annotate("text", x = t2009$dates, y = t2009$dYsim_90_cf1*1.3, label = "Q90: no MF", size = 3,color=colorspace::darken(pal2[3],0.3),hjust=0,vjust=1)+
      # geom_ribbon(aes(x=df$dates, ymax = df$dYsim_90_cf1, ymin = df$dYsim_10_cf1), alpha = 0.5,
  #             fill = pal[1], color = "transparent")+
      scale_color_identity() +
      scale_fill_identity() +
      coord_cartesian(expand = FALSE, clip = "off") +
      theme_tufte(base_size = 12, base_family = "Helvetica") +
      theme(
      axis.ticks = element_line(color = "grey50", size = 0.75),
      plot.margin=unit(c(0,1,0,0), "cm")
      )+
      labs(
      title = "Great Financial Crisis",
      subtitle = "1-year-ahead GDP growth",
      y = "Percent (%)"
      ) +
      theme(
      axis.title.x = element_blank(),
      plot.title = element_text(size = 16, face = "bold", margin=margin(0,0,0,0)),
      plot.subtitle = element_text(margin = margin(0, 0, -20, 0)),
      # plot.caption = element_text(color = "grey30")
    )+
    scale_x_date(date_breaks = "year" ,date_labels = "%Y",limits = as.Date(c("2007-01-01","2010-12-01")))+
    coord_cartesian(ylim = c(-10, 7))
p3





p4 <- p3 +
  labs(
    title = "COVID-19",
    subtitle = "1-year-ahead GDP growth",
    #caption = "Source: NBER · Graphic: Georgios Karamanis",
    y = "Percent (%)"
  ) +
  theme(
    plot.margin=unit(c(0,0,0,1), "cm")
  )+
  scale_x_date(date_breaks = "3 months" ,date_labels = "%b-%Y",limits = as.Date(c("2020-1-01","2020-10-30")))+
  annotate("text", x = t2020$dates, y = 0.5*(t2020$dYsim_10+t2020$dYsim_90), label = "Estimated", size = 3,color=colorspace::darken(pal[3],0.1),hjust=0)+
  annotate("text", x = t2020$dates, y = t2020$dYsim_10_cf1*0.70, label = "Q10: no MF", size = 3,color=colorspace::darken(pal2[5],0.1),hjust=0,vjust=1)+
  annotate("text", x = t2020$dates, y = t2020$dYsim_90_cf1*1.3, label = "Q90: no MF", size = 3,color=colorspace::darken(pal2[3],0.3),hjust=0,vjust=1)
p4


p5 <-ggplot(df) +
  geom_ribbon(aes(x=dates, ymax = dYsim_90, ymin = dYsim_10), alpha = 0.20,
              fill = pal[3], color = "transparent")+
  geom_line(data=df, aes(x=dates, y=dYsim_90_cf2), linetype='dotdash',color="steelblue4", lwd = 2.0) +
  geom_line(data=df, aes(x=dates, y=dYsim_10_cf2), linetype='dotdash',color=pal2[5], lwd = 2.0) +
  geom_line(aes(dates, dYsim_90), size = 0.5,color=colorspace::lighten(pal[3], 0.5)) +
  geom_line(aes(dates, dYsim_10), size = 0.5,color=colorspace::lighten(pal[3], 0.5)) +
  annotate("text", x = t2009$dates, y = 0.5*(t2009$dYsim_10+t2009$dYsim_90), label = "Estimated", size = 3,color=colorspace::darken(pal[3],0.1),hjust=0)+
  annotate("text", x = t2009$dates, y = t2009$dYsim_10_cf2*0.70, label = "Q10: no FF", size = 3,color=colorspace::darken(pal2[5],0.1),hjust=0,vjust=1)+
  annotate("text", x = t2009$dates, y = t2009$dYsim_90_cf2*1.3, label = "Q90: no FF", size = 3,color=colorspace::darken(pal2[3],0.3),hjust=0,vjust=1)+
  scale_color_identity() +
  scale_fill_identity() +
  coord_cartesian(expand = FALSE, clip = "off") +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(
    # legend.position = "none",
    # axis.title.x = element_blank(),
     panel.grid.minor.y = element_blank(),
     panel.grid.major.y = element_blank(),
     panel.grid.minor.x = element_blank(),
     panel.grid.major.x = element_blank(),
     axis.ticks = element_line(color = "grey50", size = 0.75)
    # axis.line.x = element_line(colour = "grey50", size = 0.5,linetype = "solid"),
     # axis.line.y = element_line(colour = "grey50", size = 0.5, linetype = "solid")
  )+
  labs(
    title = "Great Financial Crisis",
    subtitle = "1-year-ahead GDP growth",
    y = "Percent (%)"
  ) +
  theme(
    axis.title.x = element_blank(),
    plot.title = element_text(size = 16, face = "bold", margin=margin(0,0,0,0)),
    plot.subtitle = element_text(margin = margin(0, 0, -20, 0)),
    plot.caption = element_text(color = "grey30"),
    plot.margin=unit(c(0,1,0,0), "cm")
  )+
  scale_x_date(date_breaks = "year" ,date_labels = "%Y",limits = as.Date(c("2007-01-01","2010-12-01")))+
  coord_cartesian(ylim = c(-10, 7))

p5



p6 <- p5 +
  labs(
    title = "COVID-19",
    subtitle = "1-year-ahead GDP growth",
    #caption = "Source: NBER · Graphic: Georgios Karamanis",
    y = "Percent (%)"
  ) +
  theme(
    plot.margin=unit(c(0,0,0,1), "cm")
  )+
  scale_x_date(date_breaks = "3 months" ,date_labels = "%b-%Y",limits = as.Date(c("2020-1-01","2020-10-30")))+
  annotate("text", x = t2020$dates, y = 0.5*(t2020$dYsim_10+t2020$dYsim_90), label = "Estimated", size = 3,color=colorspace::darken(pal[3],0.1),hjust=0)+
  annotate("text", x = t2020$dates, y = t2020$dYsim_10_cf2*0.70, label = "Q10: no FF", size = 3,color=colorspace::darken(pal2[5],0.1),hjust=0,vjust=1)+
  annotate("text", x = t2020$dates, y = t2020$dYsim_90_cf2*1.3, label = "Q90: no FF", size = 3,color=colorspace::darken(pal2[3],0.3),hjust=0,vjust=1)

p6

p8 <- p3 + p4
p8


p9 <- p5 + p6
p9
# ggsave(paste0("quantiles", model, "_Model11_cf_2009.pdf"), p8,  width = 12, height = 8,dpi=500)
# ggsave(paste0("quantiles", model, "_Model11_cf_2020.pdf"), p9,  width = 12, height = 8,dpi=500)





bob


p1 <-ggplot(df) +
  geom_ribbon(aes(x=dates, ymax = dYsim_90, ymin = dYsim_10), alpha = 0.25,
              fill = pal[3], color = "transparent")+
  geom_line(data=df, aes(x=dates, y=dYsim_90_cf1), linetype='dashed',color=pal2[3], lwd = 2.0) +
  geom_line(data=df, aes(x=dates, y=dYsim_10_cf1), linetype='dashed',color=pal2[5], lwd = 2.0) +
  geom_line(aes(dates, dYsim_90), size = 1,color=colorspace::lighten(pal2[3], 0.2)) +
  geom_line(aes(dates, dYsim_10), size = 1,color=colorspace::lighten(pal2[5], 0.2)) +
  facet_zoom(xlim = as.Date(c("2007-01-01","2010-12-01")), horizontal = FALSE)+
  labs(
    title = "Great Financial Crisis",
    subtitle = "1-year-ahead GDP growth",
    y = "Percent (%)"
  ) +
  theme_minimal(base_family = "Source Serif Pro")+
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    plot.background = element_blank(), #element_rect(fill = "grey95", color = NA),
    strip.background = element_rect(fill = "white", color = "grey70", size = 1),
    panel.background = element_blank(), #element_rect(fill = "grey90", color = NA),
    axis.title = element_blank(),
    axis.text.x = element_text(face = "bold"),
    plot.margin = margin(20, 20, 20, 20),
    # plot.title = element_text(hjust = 0.5, size = 23, family = "Fira Sans", face = "bold", color = "#105585"),
    # plot.subtitle = element_text(hjust = 0.5, size = 16, family = "Fira Sans", face = "bold", color = "#105585"),
    plot.caption = element_text(family = "Fira Sans", color = "grey40"),
    panel.grid.major.y = element_blank(), #element_line(color = "grey50", size = 0.125),
    panel.grid.minor.y = element_blank(), #element_line(color = "grey50", size = 0.05),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(size = 12)
  )

p1


p2 <-ggplot(df) +
  geom_ribbon(aes(x=dates, ymax = dYsim_90, ymin = dYsim_10), alpha = 0.25,
              fill = pal[3], color = "transparent")+
  geom_line(data=df, aes(x=dates, y=dYsim_90_cf1), linetype='dashed',color=pal2[3], lwd = 2.0) +
  geom_line(data=df, aes(x=dates, y=dYsim_10_cf1), linetype='dashed',color=pal2[5], lwd = 2.0) +
  geom_line(aes(dates, dYsim_90), size = 1,color=colorspace::lighten(pal2[3], 0.2)) +
  geom_line(aes(dates, dYsim_10), size = 1,color=colorspace::lighten(pal2[5], 0.2)) +
  facet_zoom(xlim = as.Date(c("2020-01-01","2020-09-30")), horizontal = FALSE)+
  labs(
    title = "COVID-19",
    subtitle = "1-year-ahead GDP growth",
    y = "Percent (%)"
  ) +
  theme_minimal(base_family = "Source Serif Pro")+
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    plot.background = element_blank(), #element_rect(fill = "grey95", color = NA),
    strip.background = element_rect(fill = "white", color = "grey70", size = 1),
    panel.background = element_blank(), #element_rect(fill = "grey90", color = NA),
    axis.title = element_blank(),
    axis.text.x = element_text(face = "bold"),
    plot.margin = margin(20, 20, 20, 20),
    # plot.title = element_text(hjust = 0.5, size = 23, family = "Fira Sans", face = "bold", color = "#105585"),
    # plot.subtitle = element_text(hjust = 0.5, size = 16, family = "Fira Sans", face = "bold", color = "#105585"),
    plot.caption = element_text(family = "Fira Sans", color = "grey40"),
    panel.grid.major.y = element_blank(), #element_line(color = "grey50", size = 0.125),
    panel.grid.minor.y = element_blank(), #element_line(color = "grey50", size = 0.05),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(size = 12)
  )

p2








ggplot(df) +
  geom_line(aes(dates, dYsim_10), size = 1) +
  facet_zoom(xlim = as.Date(c("2007-01-01","2010-12-01")), horizontal = FALSE)
  labs(
    title = "Real Digital Infrastructure Investment",
    subtitle = "In chained (2012) dollars",
    caption = "Source: Bureau of Economic Analysis · Graphic: Georgios Karamanis"
  ) +
  theme_minimal(base_family = "Source Serif Pro") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    plot.background = element_rect(fill = "grey95", color = NA),
    strip.background = element_rect(fill = "grey90", color = "grey70", size = 1),
    panel.background = element_rect(fill = "grey90", color = NA),
    axis.title = element_blank(),
    axis.text.x = element_text(face = "bold"),
    plot.margin = margin(20, 20, 20, 20),
    plot.title = element_text(hjust = 0.5, size = 23, family = "Fira Sans", face = "bold", color = "#105585"),
    plot.subtitle = element_text(hjust = 0.5, size = 16, family = "Fira Sans", face = "bold", color = "#105585"),
    plot.caption = element_text(family = "Fira Sans", color = "grey40"),
    panel.grid.major.y = element_line(color = "grey50", size = 0.125),
    panel.grid.minor.y = element_blank(), #element_line(color = "grey50", size = 0.05),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(size = 12)
  )

