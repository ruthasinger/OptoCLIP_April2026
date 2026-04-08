#load libraries 
library(tidyverse)
library(readABF)
library(zoo)

##############################################################################################################
#setup the directories
Base_directory=file.path("~/ruthasinger_github/OptoCLIP_April2026/Figure2_FigureS2_Ephys")
list.files(Base_directory)

Data_directory=file.path(Base_directory,"Data")
list.files(Data_directory)

Outdirectory=file.path(Base_directory,paste("Output",Sys.Date(),sep="_"))
dir.create(Outdirectory,showWarnings = TRUE)

##############################################################################################################
#current clamp data showing health of neurons that we are doing opto on

`%||%` <- function(a, b) if (!is.null(a)) a else b

read_cc_vm_only_abf <- function(abf_path, vm_pattern = "\\[mV\\]") {
  abf <- readABF(abf_path)
  nsweeps <- length(abf$data)
  
  dfs <- lapply(seq_len(nsweeps), function(s) {
    df <- as.data.frame(abf, sweep = s) |> as.data.frame()
    nm <- names(df)
    
    time_col <- grep("Time", nm, value = TRUE)[1]
    vm_col   <- grep(vm_pattern, nm, value = TRUE)[1]
    
    if (any(is.na(c(time_col, vm_col)))) {
      stop("Could not find Time/mV in ", basename(abf_path),
           "\nColumns are: ", paste(nm, collapse = ", "))
    }
    
    tibble(Time = df[[time_col]], sweep = s, mV = df[[vm_col]])
  })
  
  bind_rows(dfs)
}

# --- Extract pA + stim window from donor file's AO #0 [pA] ---
extract_sweep_pA_map_from_AO <- function(abf_path,
                                         ao_col = "AO #0 [pA]",
                                         min_step_pA = 5) {
  abf <- readABF(abf_path)
  nsweeps <- length(abf$data)
  
  lapply(seq_len(nsweeps), function(s) {
    df <- as.data.frame(abf, sweep = s) |> as.data.frame()
    nm <- names(df)
    
    time_col <- grep("Time", nm, value = TRUE)[1]
    if (is.na(time_col) || !(ao_col %in% nm)) {
      stop("Could not find Time or AO column in ", basename(abf_path),
           "\nColumns are: ", paste(nm, collapse = ", "))
    }
    
    t  <- df[[time_col]]
    ao <- df[[ao_col]]
    
    # baseline from first 20% of sweep
    n <- length(ao)
    base_idx <- seq_len(max(10, floor(0.2 * n)))
    base <- median(ao[base_idx], na.rm = TRUE)
    
    idx <- which(abs(ao - base) > min_step_pA)
    
    if (length(idx) == 0) {
      tibble(sweep = s, pA = 0, stim_start = NA_real_, stim_end = NA_real_)
    } else {
      tibble(
        sweep = s,
        pA = round(median(ao[idx], na.rm = TRUE), 1),
        stim_start = t[min(idx)],
        stim_end   = t[max(idx)]
      )
    }
  }) %>% bind_rows()
}

# --- Option A spike counting: exclude NA windows; join back so 0 pA gets 0 spikes ---
count_spikes_in_window <- function(cc_df, stim_df,
                                   spike_thresh = 0,
                                   refractory = 0.001) {
  
  stim_df_valid <- stim_df %>% filter(!is.na(stim_start), !is.na(stim_end))
  
  spikes <- cc_df %>%
    inner_join(stim_df_valid, by = "sweep") %>%
    filter(Time >= stim_start, Time <= stim_end) %>%
    group_by(sweep) %>%
    arrange(Time) %>%
    mutate(
      above = mV > spike_thresh,
      spike = above & !lag(above, default = FALSE)
    ) %>%
    filter(spike) %>%
    mutate(
      isi = Time - lag(Time),
      valid = is.na(isi) | isi > refractory
    ) %>%
    filter(valid) %>%
    summarise(n_APs = n(), .groups = "drop")
  
  FI_df <- stim_df %>%
    left_join(spikes, by = "sweep") %>%
    mutate(n_APs = replace_na(n_APs, 0L)) %>%
    arrange(pA)
  
  FI_df
}

# --- Per-cell analysis wrapper (works for Vm-only + donor file) ---
analyze_one_cc_cell <- function(abf_path, sweep_pA_map,
                                spike_thresh = 0, refractory = 0.001) {
  
  cc_df <- read_cc_vm_only_abf(abf_path) %>%
    left_join(sweep_pA_map %>% dplyr::select(sweep, pA), by = "sweep")
  
  FI_df <- count_spikes_in_window(
    cc_df = cc_df,
    stim_df = sweep_pA_map,
    spike_thresh = spike_thresh,
    refractory = refractory
  ) %>%
    mutate(cell = basename(abf_path))
  
  FI_df
}

get_rheobase_for_cell <- function(FI_df) {
  rb <- FI_df %>%
    filter(n_APs >= 1) %>%
    summarise(rb = min(pA)) %>%
    pull(rb)
  if (length(rb) == 0 || !is.finite(rb)) NA_real_ else rb
}

donor_abf <- file.path(Data_directory, "CurrentClamp",
                       "08112021_Camk2a_ChR2_Cell1_cc200pA.abf")

other_vm_only <- c(
  "optoCamK2a_Aug11_cc200pAsweep.abf",
  "optoCamK2a_Aug11_cc200pAsweep_postLED.abf"
)

other_vm_only <- file.path(Data_directory, "CurrentClamp", other_vm_only)

# Include donor file as a data file too:
all_cells <- c(other_vm_only, donor_abf)

#Build protocol map (pA + window) from donor file
sweep_pA_map <- extract_sweep_pA_map_from_AO(donor_abf, ao_col = "AO #0 [pA]", min_step_pA = 5)

#Get FI per cell
FI_all_cells <- lapply(all_cells, analyze_one_cc_cell, sweep_pA_map = sweep_pA_map) %>%
  bind_rows()

#Average FI across cells
FI_summary <- FI_all_cells %>%
  group_by(pA) %>%
  summarise(
    mean_APs = mean(n_APs),
    sem_APs  = sd(n_APs) / sqrt(n()),
    n_cells  = n(),
    .groups = "drop"
  )

Rheostat_mean_plot <- ggplot(FI_summary, aes(x = pA, y = mean_APs)) +
  # SEM shading
  geom_ribbon(
    aes(ymin = mean_APs - sem_APs,
        ymax = mean_APs + sem_APs),
    fill = "gray70",
    alpha = 0.4
  ) +
  
  # mean line
  geom_line(linewidth = 0.6) +
  
  # mean points
  geom_point(size = 1.5) +
  
  # rheobase line
  geom_vline(
    xintercept = 50,
    linetype = "dashed",
    colour = "green",
    linewidth = 0.5
  ) +
  
  theme_classic(base_size = 9) +
  labs(
    x = "Injected current (pA)",
    y = "Mean # APs ± SEM"
  )
Rheostat_mean_plot
ggsave(file.path(Outdirectory, "Rheostat_mean_plot_across_3_cells.pdf"), Rheostat_mean_plot, device = "pdf", width = 3, height = 2, units = "in")

get_rheobase_for_cell <- function(FI_df) {
  rb <- FI_df %>%
    filter(n_APs >= 1) %>%
    summarise(rb = min(pA)) %>%
    pull(rb)
  if (length(rb) == 0 || !is.finite(rb)) NA_real_ else rb
}

make_rheobase_color_map <- function(sweep_pA_map, rheobase_pA) {
  
  sweep_order <- sweep_pA_map %>% arrange(pA) %>% pull(sweep)
  
  cols <- setNames(
    gray.colors(length(sweep_order), start = 0.1, end = 0.7),
    as.character(sweep_order)
  )
  
  if (!is.na(rheobase_pA)) {
    rb_sweep <- sweep_pA_map %>%
      filter(pA == rheobase_pA) %>%
      pull(sweep)
    
    if (length(rb_sweep) > 0) {
      cols[as.character(rb_sweep[1])] <- "green"
    }
  }
  
  cols
}

plot_Vm_one_cell <- function(abf_path, sweep_pA_map,
                             xlim = c(5,6),
                             spike_thresh = 0,
                             refractory = 0.001) {
  
  cc_df <- read_cc_vm_only_abf(abf_path) %>%
    left_join(sweep_pA_map %>% dplyr::select(sweep, pA), by = "sweep")
  
  FI_df <- count_spikes_in_window(
    cc_df = cc_df,
    stim_df = sweep_pA_map,
    spike_thresh = spike_thresh,
    refractory = refractory
  )
  
  rheobase_pA <- get_rheobase_for_cell(FI_df)
  
  sweep_order <- sweep_pA_map %>% arrange(pA) %>% pull(sweep)
  cc_df <- cc_df %>% mutate(sweep = factor(sweep, levels = sweep_order))
  
  cols <- make_rheobase_color_map(sweep_pA_map, rheobase_pA)
  
  Vm_plot <- ggplot(cc_df, aes(x = Time, y = mV, color = sweep, group = sweep)) +
    coord_cartesian(xlim = xlim) +
    geom_line(linewidth = 1, show.legend = FALSE) +
    theme_void() +
    scale_color_manual(values = cols, drop = FALSE) +
    annotate("segment", x = 5.9, y= -80, xend = 6.0, yend= -80, color="black", linewidth=1) +
    annotate("segment", x = 6.0, y= -80, xend = 6.0, yend= -70, color="black", linewidth=1) +
    ggtitle(paste0(
      basename(abf_path),
      ifelse(is.na(rheobase_pA),
             " (no spikes)",
             paste0(" | rheobase = ", rheobase_pA, " pA"))
    ))
  print(Vm_plot)
  ggsave(
    filename = file.path(Outdirectory, "Vm_plot_08112021_Camk2a_ChR2_Cell1_cc200pA_scale_100ms_10mV.pdf"),
    plot = Vm_plot, device = "pdf", width = 5, height = 4, units = "in"
  )
  
  invisible(list(plot = Vm_plot, rheobase_pA = rheobase_pA, cols = cols))
}

plot_pA_one_cell <- function(abf_path,
                             sweep_pA_map,
                             rheobase_pA,
                             ao_col = "AO #0 [pA]",
                             xlim = c(5,6)) {
  
  abf <- readABF(abf_path)
  nsweeps <- length(abf$data)
  
  df <- lapply(seq_len(nsweeps), function(s) {
    d <- as.data.frame(abf, sweep = s)
    tibble(
      Time  = d[["Time [s]"]],
      sweep = s,
      pA    = d[[ao_col]]
    )
  }) %>% bind_rows()
  
  sweep_order <- sweep_pA_map %>% arrange(pA) %>% pull(sweep)
  df <- df %>% mutate(sweep = factor(sweep, levels = sweep_order))
  
  # SAME exact color scheme as Vm plot
  cols <- make_rheobase_color_map(sweep_pA_map, rheobase_pA)
  
  pA_plot <- ggplot(df, aes(x = Time, y = pA, group = sweep, color = sweep)) +
    coord_cartesian(xlim = xlim) +
    geom_line(linewidth = 0.8, show.legend = FALSE) +
    theme_void(base_size = 9) +
    scale_color_manual(values = cols, drop = FALSE) +
    labs(
      title = paste0(basename(abf_path), " — injected current"),
      x = "Time (s)",
      y = "Injected current (pA)"
    ) +
    annotate("segment", x = 5.9, y= -50, xend = 6.0, yend= -50, color="black", linewidth=1) +
    annotate("segment", x = 6.0, y= -50, xend = 6.0, yend= -25, color="black", linewidth=1)
  
  print(pA_plot)
  ggsave(
    filename = file.path(Outdirectory, "pA_plot_08112021_Camk2a_ChR2_Cell1_cc200pA_scale_100ms_25pA.pdf"),
    plot = pA_plot, device = "pdf", width = 5, height = 4, units = "in"
  )
  
  invisible(pA_plot)
}

# usage
FI_df_donor <- analyze_one_cc_cell(donor_abf, sweep_pA_map)
rheobase_pA_donor <- get_rheobase_for_cell(FI_df_donor)

Vm_res <- plot_Vm_one_cell(all_cells[3], sweep_pA_map)
pA_plot <- plot_pA_one_cell(all_cells[3], sweep_pA_map, rheobase_pA_donor)

##############################################################################################################
#showing spiking at different Hz LED, 1 Hz, 5Hz, and 10 Hz

read_abf_trace <- function(abf_file, value_name, channel = 1, gain_divide = 1, n_max = NULL) {
  x <- readABF(abf_file)
  
  mat <- x$data[[1]]
  
  dt <- x$samplingIntervalInSec * ncol(mat)
  
  df <- data.frame(
    Time = seq(0, by = dt, length.out = nrow(mat)),
    value = mat[, channel] / gain_divide
  )
  
  if (!is.null(n_max)) {
    df <- df[seq_len(min(n_max, nrow(df))), , drop = FALSE]
  }
  
  colnames(df)[2] <- value_name
  df
}

# Control files
LED_1Hz_Control <- read_abf_trace(
  file.path(Data_directory, "DifferentHz", "020222_Patch_Control_1Hz_bluelight.abf"),
  "LED_1Hz_Control"
)

LED_5Hz_Control <- read_abf_trace(
  file.path(Data_directory, "DifferentHz", "020222_Patch_Control_5Hz_bluelight.abf"),
  "LED_5Hz_Control"
)

LED_10Hz_Control <- read_abf_trace(
  file.path(Data_directory, "DifferentHz", "020222_Patch_Control_10Hz_bluelight.abf"),
  "LED_10Hz_Control"
)

# ChR2 files
LED_1Hz_ChR2 <- read_abf_trace(
  file.path(Data_directory, "DifferentHz", "033123_Patch_ChR2_1Hz_bluelight.abf"),
  "LED_1Hz_ChR2",
  channel = 1,
  gain_divide = 20,
  n_max = 210000
)

LED_5Hz_ChR2 <- read_abf_trace(
  file.path(Data_directory, "DifferentHz", "033123_Patch_ChR2_5Hz_bluelight.abf"),
  "LED_5Hz_ChR2",
  channel = 1,
  gain_divide = 20,
  n_max = 210000
)

LED_10Hz_ChR2 <- read_abf_trace(
  file.path(Data_directory, "DifferentHz", "033123_Patch_ChR2_10Hz_bluelight.abf"),
  "LED_10Hz_ChR2",
  channel = 1,
  gain_divide = 20,
  n_max = 210000
)

#1Hz
segment_data1hz = data.frame(
  x = seq(1.25,20,2),
  xend=seq(1.25,20,2),
  y = -75,
  yend = -72)
nrow(segment_data1hz)

OneHz_plot=ggplot() +
  geom_line(data = LED_1Hz_ChR2, aes(x = Time, y = LED_1Hz_ChR2, color = "ChR2"), linewidth=.5) +
  geom_line(data = LED_1Hz_Control, aes(x = Time, y = LED_1Hz_Control, color = "Control"), linewidth=.5) +
  scale_color_manual(name = "AAV",
                     values = c("ChR2" = "red", "Control" = "black")) +
  geom_smooth(method = "loess") +
  xlim(0,20) +
  geom_segment(data = segment_data1hz, aes(x = x, y = y, xend = xend, yend = yend), color="steelblue2",linewidth=.5) +
  annotate("segment",x = 15, y =40, xend = 17, yend =40, color="black",linewidth=.5) + 
  annotate("segment", x = 17, y =40, xend =17, yend =50, color="black",linewidth=.5) +
  theme_void() +
  theme(legend.position="none")
OneHz_plot
ggsave(filename = file.path(Outdirectory,"LED_1Hz_ChR2andControl_scale_1s_10mV.pdf"), plot= OneHz_plot, device='pdf', width=1.7,height=2.25,unit="in")

OneHz_plot_zoom=ggplot() +
  geom_line(data = LED_1Hz_ChR2, aes(x = Time, y = LED_1Hz_ChR2, color = "ChR2"), linewidth=.5) +
  geom_line(data = LED_1Hz_Control, aes(x = Time, y = LED_1Hz_Control, color = "Control"), linewidth=.5) +
  scale_color_manual(name = "AAV",
                     values = c("ChR2" = "red", "Control" = "black")) +
  geom_smooth(method = "loess") +
  xlim(4,6) +
  annotate("segment",x = 5.7, y =40, xend = 5.9, yend =40, color="black",linewidth=.5) + 
  annotate("segment", x = 5.9, y =40, xend =5.9, yend =50, color="black",linewidth=.5) +
  geom_segment(data = segment_data1hz, aes(x = x, y = y, xend = xend, yend = yend), color="steelblue2",linewidth=.5) +
  theme_void() +
  theme(legend.position="none")
OneHz_plot_zoom
ggsave(filename = file.path(Outdirectory,"LED_1Hz_ChR2andControl_scale_100ms_10mV_zoom_4to6.pdf"), plot= OneHz_plot_zoom, device='pdf', width=1.7,height=2.25,unit="in")

#5hz
segment_data5hz = data.frame(
  x = seq(1.25,20,.4),
  xend=seq(1.25,20,.4),
  y = -75,
  yend = -72)

FiveHz_plot=ggplot() +
  geom_line(data = LED_5Hz_ChR2, aes(x = Time, y = LED_5Hz_ChR2, color = "ChR2"), linewidth=.5) +
  geom_line(data = LED_5Hz_Control, aes(x = Time, y = LED_5Hz_Control, color = "Control"), linewidth=.5) +
  scale_color_manual(name = "AAV",
                     values = c("ChR2" = "red", "Control" = "black")) +
  geom_smooth(method = "loess") +
  xlim(0,20) +
  geom_segment(data = segment_data5hz, aes(x = x, y = y, xend = xend, yend = yend), color="steelblue2",linewidth=.5) +
  annotate("segment",x = 15, y =40, xend = 17, yend =40, color="black",linewidth=.5) + 
  annotate("segment", x = 17, y =40, xend =17, yend =50, color="black",linewidth=.5) +
  theme_void() +
  theme(legend.position="none")
FiveHz_plot
ggsave(filename = file.path(Outdirectory,"LED_5Hz_ChR2andControl_scale_1s_10mV.pdf"), plot= FiveHz_plot, device='pdf',width=1.7,height=2.25,unit="in")

FiveHz_plot_zoom=ggplot() +
  geom_line(data = LED_5Hz_ChR2, aes(x = Time, y = LED_5Hz_ChR2, color = "ChR2"), linewidth=.5) +
  geom_line(data = LED_5Hz_Control, aes(x = Time, y = LED_5Hz_Control, color = "Control"), linewidth=.5) +
  scale_color_manual(name = "AAV",
                     values = c("ChR2" = "red", "Control" = "black")) +
  geom_smooth(method = "loess") +
  xlim(4,6) +
  annotate("segment",x = 5.7, y =40, xend = 5.9, yend =40, color="black",linewidth=.5) + 
  annotate("segment", x = 5.9, y =40, xend =5.9, yend =50, color="black",linewidth=.5) +
  geom_segment(data = segment_data5hz, aes(x = x, y = y, xend = xend, yend = yend), color="steelblue2",linewidth=.5) +
  theme_void() +
  theme(legend.position="none")
FiveHz_plot_zoom
ggsave(filename = file.path(Outdirectory,"LED_5Hz_ChR2andControl_scale_100ms_10mV_zoom_4to6.pdf"), plot= FiveHz_plot_zoom, device='pdf',width=1.7,height=2.25,unit="in")

#10Hz
segment_data10hz = data.frame(
  x = seq(1.25,20,.2),
  xend=seq(1.25,20,.2),
  y = -75,
  yend = -72)

TenHz_plot=ggplot() +
  geom_line(data = LED_10Hz_ChR2, aes(x = Time, y = LED_10Hz_ChR2, color = "ChR2"), linewidth=.5) +
  geom_line(data = LED_10Hz_Control, aes(x = Time, y = LED_10Hz_Control, color = "Control"), linewidth=.5) +
  scale_color_manual(name = "AAV",
                     values = c("ChR2" = "red", "Control" = "black")) +
  geom_smooth(method = "loess") +
  xlim(0,20) +
  geom_segment(data = segment_data10hz, aes(x = x, y = y, xend = xend, yend = yend), color="steelblue2",linewidth=.5) +
  annotate("segment",x = 15, y =40, xend = 17, yend =40, color="black",linewidth=.5) + 
  annotate("segment", x = 17, y =40, xend =17, yend =50, color="black",linewidth=.5) +
  theme_void() +
  theme(legend.position="none")
TenHz_plot
ggsave(filename = file.path(Outdirectory,"LED_10Hz_ChR2andControl_scale_1s_10mV.pdf"), plot= TenHz_plot, device='pdf', width=1.7,height=2.25,unit="in")

TenHz_plot_zoom=ggplot() +
  geom_line(data = LED_10Hz_ChR2, aes(x = Time, y = LED_10Hz_ChR2, color = "ChR2"), linewidth=.5) +
  geom_line(data = LED_10Hz_Control, aes(x = Time, y = LED_10Hz_Control, color = "Control"), linewidth=.5) +
  scale_color_manual(name = "AAV",
                     values = c("ChR2" = "red", "Control" = "black")) +
  geom_smooth(method = "loess") +
  xlim(4,6) +
  annotate("segment",x = 5.7, y =40, xend = 5.9, yend =40, color="black",linewidth=.5) + 
  annotate("segment", x = 5.9, y =40, xend =5.9, yend =50, color="black",linewidth=.5) +
  geom_segment(data = segment_data10hz, aes(x = x, y = y, xend = xend, yend = yend), color="steelblue2",linewidth=.5) +
  theme_void() +
  theme(legend.position="none")
TenHz_plot_zoom
ggsave(filename = file.path(Outdirectory,"LED_10Hz_ChR2andControl_scale_100ms_10mV_zoom_4to6.pdf"), plot= TenHz_plot_zoom, device='pdf', width=1.7,height=2.25,unit="in")

##############################################################################################################
#showing spiking with Opto-CLIP paradigm- 5 trains of 30" light on (5 Hz pulses), 30" light off

abf_files <- c(
  "042325_ChR2_Neuron_1.abf",
  "042325_ChR2_Neuron_2.abf",
  "042325_ChR2_Neuron_3.abf",
  "042325_Control_Neuron_1.abf",
  "042325_Control_Neuron_2.abf"
)

labels <- sub("\\.abf$", "", abf_files)

read_all_sweeps <- function(filename) {
  abf <- readABF(file.path(Data_directory, "5trains5Hz",filename))
  nsweeps <- length(abf$data)
  
  sweep_dfs <- lapply(seq_len(nsweeps), function(s) {
    df <- as.data.frame(abf, sweep = s)
    colnames(df) <- c("Time", "Signal")
    df$Sweep <- s
    return(df)
  })
  
  do.call(rbind, sweep_dfs)
}

Neuron_Data <- setNames(lapply(abf_files, read_all_sweeps), labels)

Neuron_Data <- Map(function(df, name) {
  df$Condition <- if (grepl("Control", name)) "Control" else "ChR2"
  df$Neuron <- name
  df
}, Neuron_Data, labels)

all_neurons_df <- do.call(rbind, Neuron_Data)

neurons_to_plot <- c("042325_ChR2_Neuron_3",
                     "042325_Control_Neuron_2") 

plot_df <- all_neurons_df %>%
  filter(Neuron %in% neurons_to_plot, Sweep %in% 1:5) %>%
  group_by(Neuron, Sweep) %>%
  mutate(Signal_Smooth = zoo::rollmean(Signal, 50, fill = NA)) %>%
  ungroup()

stim_start_df <- plot_df %>%
  group_by(Neuron, Sweep) %>%
  filter(Signal > 20) %>%
  summarise(first_spike_time = min(Time, na.rm = TRUE), .groups = "drop") %>%
  mutate(stim_start = first_spike_time)

pulse_df <- stim_start_df %>%
  rowwise() %>%
  mutate(pulses = list(seq(stim_start, stim_start + 29.8, by = 0.2))) %>%
  unnest(cols = c(pulses)) %>%
  rename(Time = pulses)

signal_min <- min(plot_df$Signal, na.rm = TRUE)
tick_bottom <- signal_min - 0.05 * abs(signal_min)
tick_top <- tick_bottom + 0.05 * abs(signal_min)

tick_df <- pulse_df %>%
  mutate(x = Time, xend = Time, y = tick_bottom, yend = tick_top)

plot_df <- plot_df %>%
  left_join(stim_start_df, by = c("Neuron", "Sweep"))

scale_df <- stim_start_df %>%
  mutate(
    scale_x = stim_start - 2.2,
    scale_xend = stim_start - 1.2,
    scale_y = tick_bottom + 0.1 * abs(tick_bottom),
    scale_yend = scale_y + 10
  )

FiveHz_FiveTrain_plot <- ggplot(
  plot_df,
  aes(x = Time, y = Signal,
      group = interaction(Neuron, Sweep),
      color = Condition)
) +
  geom_line(linewidth = 0.25, alpha = 0.9) +
  
  geom_segment(
    data = tick_df,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE,
    color = "dodgerblue",
    linewidth = 0.25
  ) +
  
  geom_segment(
    data = scale_df,
    aes(x = scale_x, xend = scale_xend,
        y = scale_y, yend = scale_y),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.35
  ) +
  
  geom_segment(
    data = scale_df,
    aes(x = scale_x, xend = scale_x,
        y = scale_y, yend = scale_yend),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.35
  ) +
  
  scale_color_manual(values = c(
    "ChR2" = "red3",
    "Control" = "black"
  )) +
  
  facet_wrap(
    ~Sweep,
    ncol = 1,
    strip.position = "left"
  ) +
  
  theme_void(base_size = 8) +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank()
  )
FiveHz_FiveTrain_plot
ggsave(filename = file.path(Outdirectory,"FiveHz_FiveTrain_plot_scale_1s_10mV.png"), plot= FiveHz_FiveTrain_plot, device='png', dpi = 600, width=7, height=3.5)

###########################
#making summary plot
abf_files <- c(
  "042325_ChR2_Neuron_1.abf",
  "042325_ChR2_Neuron_2.abf",
  "042325_ChR2_Neuron_3.abf",
  "042325_Control_Neuron_1.abf",
  "042325_Control_Neuron_2.abf"
)

labels <- sub("\\.abf$", "", abf_files)

read_all_sweeps <- function(filename) {
  abf <- readABF(file.path(Data_directory, "5trains5Hz",filename))
  nsweeps <- length(abf$data)
  
  bind_rows(lapply(seq_len(nsweeps), function(s) {
    df <- as.data.frame(abf, sweep = s)
    colnames(df) <- c("Time", "Signal")  # Time [s], IN 0 [mV]
    df$Sweep <- s
    df
  }))
}

Neuron_Data <- setNames(lapply(abf_files, read_all_sweeps), labels)

Neuron_Data <- Map(function(df, name) {
  df$Condition <- if (grepl("Control", name)) "Control" else "ChR2"
  df$Neuron <- name
  df
}, Neuron_Data, labels)

all_neurons_df <- bind_rows(Neuron_Data)

stim_start_df_all <- all_neurons_df %>%
  filter(Sweep %in% 1:5) %>%
  group_by(Neuron, Sweep) %>%
  summarise(
    stim_start = min(Time, na.rm = TRUE),   # usually 0
    stim_end   = min(Time, na.rm = TRUE) + 30,
    .groups = "drop"
  ) %>%
  left_join(all_neurons_df %>% distinct(Neuron, Condition), by = "Neuron")

detect_spikes_times <- function(Time, Signal, thresh = 0, refrac = 0.003) {
  idx <- which(Signal[-1] >= thresh & Signal[-length(Signal)] < thresh) + 1
  t <- Time[idx]
  if (length(t) > 1) t <- t[c(TRUE, diff(t) >= refrac)]
  t
}

spike_thresh <- -10    # 0 mV is usually safer than 20 mV
refrac_s     <- 0.05

fr_df <- all_neurons_df %>%
  filter(Sweep %in% 1:5) %>%
  inner_join(stim_start_df_all, by = c("Neuron", "Sweep", "Condition")) %>%
  group_by(Neuron, Condition, Sweep) %>%
  group_modify(~{
    st  <- unique(.x$stim_start)
    end <- unique(.x$stim_end)
    
    spk <- detect_spikes_times(.x$Time, .x$Signal, thresh = spike_thresh, refrac = refrac_s)
    
    n_light <- sum(spk >= st  & spk < end)
    n_off   <- sum(spk >= end & spk < (end + 30))
    
    tibble(
      n_spikes_light = n_light,
      FR_light = n_light / 30,
      n_spikes_off   = n_off,
      FR_off   = n_off / 30
    )
  }) %>%
  ungroup()

fr_df %>% arrange(Condition, Neuron, Sweep)

fr_df2 <- fr_df %>%
  mutate(
    Sweep = factor(Sweep, levels = 1:5),
    Neuron = factor(Neuron),
    Condition = factor(Condition, levels = c("Control", "ChR2"))
  )

fr_long <- fr_df2 %>%
  tidyr::pivot_longer(c(FR_light, FR_off), names_to = "Epoch", values_to = "FR_Hz") %>%
  mutate(Epoch = recode(Epoch, FR_light = "Light (0–30s)", FR_off = "Post (30–60s)"))

fr_long <- fr_df %>%
  mutate(
    Sweep = factor(Sweep, levels = 1:5),
    Condition = factor(Condition, levels = c("ChR2", "Control")),
    Neuron = factor(Neuron)
  ) %>%
  tidyr::pivot_longer(
    cols = c(FR_light, FR_off),
    names_to = "Epoch",
    values_to = "FR_Hz"
  ) %>%
  mutate(
    Epoch = factor(
      Epoch,
      levels = c("FR_light", "FR_off"),
      labels = c("Light ON (0–30 s)", "Light OFF (30–60 s)")
    )
  )

fr_long <- fr_df %>%
  mutate(
    Sweep = factor(Sweep, levels = 1:5),
    Condition = factor(Condition, levels = c("ChR2", "Control")),
    Neuron = factor(Neuron)
  ) %>%
  pivot_longer(
    cols = c(FR_light, FR_off),
    names_to = "Epoch",
    values_to = "FR_Hz"
  ) %>%
  mutate(
    Epoch = factor(
      Epoch,
      levels = c("FR_light", "FR_off"),
      labels = c("Light ON (0–30 s)", "Light OFF (30–60 s)")
    )
  )

bg_df <- fr_long %>%
  distinct(Epoch, Sweep) %>%
  filter(Epoch == "Light ON (0–30 s)") %>%
  mutate(
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = Inf
  )

mean_df <- fr_long %>%
  group_by(Sweep, Epoch, Condition) %>%
  summarise(mean_FR = mean(FR_Hz, na.rm = TRUE), .groups = "drop")

mean_df <- mean_df %>%
  mutate(
    x_num = as.numeric(Condition),
    x    = x_num - 0.10,
    xend = x_num + 0.10
  )

p_side_by_side_bar <- ggplot(fr_long, aes(x = Condition, y = FR_Hz, fill = Condition)) +
  
  geom_rect(
    data = bg_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "lightblue",
    alpha = 0.25
  ) +
  
  stat_summary(
    fun = mean,
    geom = "col",
    width = 0.55,
    alpha = 0.6,
    color = "black"
  ) +
  
  geom_point(
    position = position_jitter(width = 0.08),
    size = 2.2,
    color = "black"
  ) +
  
  facet_grid(Sweep ~ Epoch) +
  scale_fill_manual(values = c("ChR2" = "red3", "Control" = "black")) +
  
  scale_y_continuous(
    breaks = c(0, 2.5, 5),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  
  coord_cartesian(ylim = c(-0.05, 5)) +
  
  theme_classic(base_size = 13) +
  labs(x = NULL, y = "Firing rate (Hz)") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

p_side_by_side_bar
ggsave(filename = file.path(Outdirectory,"FiveHz_FiveTrain_Hz_Allneurons.pdf"), plot= p_side_by_side_bar, device='pdf', width=5, height=7)
