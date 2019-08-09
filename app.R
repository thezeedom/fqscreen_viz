# CONTROL WIDGETS
# References:
# https://github.com/aagarw30/R-Shinyapp-Tutorial/tree/master/fileinput
# https://www.youtube.com/watch?v=HPZSunrSo5M

#sink("log_console.txt")

.packages = c("shiny",
              "shinythemes",
              "plyr",
              "dplyr",
              "reshape2",
              "DT",
              "ggplot2")

sapply(.packages, require, character.only = TRUE)
source("fqscreen_functions_MASTER.R")
#source("fqscreen_functions_SCRATCH.R")
options(scipen = 50)

#
#
#
# ------------ SHINY UI ----------- #
ui = fluidPage(
  
  theme = shinytheme("cerulean"),
  
  tags$head(
    tags$style(HTML("#next_plt { float: right; margin-left: 6px; color: #CB3929; border-color: #CB3929; }")),
    tags$style(HTML("#last_plt { float: right; color: #CB3929; border-color: #CB3929; }"))
  ),
  
  titlePanel("FastQ Screen Viz"),
  sidebarLayout(
    # Sidebar panel
    sidebarPanel(
      width = 3,
      fileInput("files", label = "Choose files:",
                multiple = TRUE,
                accept = c("text", ".txt")),
      helpText("Required file format: FastQ Screen output, .txt"),
      helpText("To select multiple files, please run this app in browser.")
    ),
    # Main Panel
    mainPanel(
      uiOutput("tb")
    )
  )
)

#
#
#
# ---------- SHINY SERVER ---------- #
server = function(input, output) {
  
  # Get data from fileInput, ID = "files"
  data = reactive({
    all_files = input$files
    if(is.null(all_files)){return()}
    # Function to create plotting df from all files
    dfs = df_from_all_files(all_files)
    df = dfs[[1]]
    return(df)
  })
  
  num_data = reactive({
    all_files = input$files
    if(is.null(all_files)){return()}
    # Function to create plotting df from all files
    dfs = df_from_all_files(all_files)
    df = dfs[[2]]
    return(df)
  })
  
  # Reactive output
  # Summary of uploaded data
  output$files_df = renderTable({
    if(is.null(data())){return()}
    input$files
  })
  
  # DF data for plotting
  output$data_summary = DT::renderDataTable({
    if(is.null(data())){return()}
    df = data()
  })
  
  # Download DF
  output$data_dl = downloadHandler(
    filename = function(){"fqscreen_viz_data.csv"}, 
    content = function(fname){
      write.csv(data(), fname, row.names = FALSE)
    }
  )
  
  samples = reactive({
    unique(data()$sample_ord)
  })
  
  is_single = reactive({
    all(data()$read=="Single-End")
  })
  
  # Select sample to plot
  observeEvent(input$which_sample, {
    if(is.null(data())){return()}
    cat("------ SAMPLE SELECTED ------", sep = "\n")
    samp_len <<- length(samples())
    cat("Single-End reads is:", is_single(), "\n")
    cat("Number of samples is: ", samp_len, "\n")
    samp_name <<- input$which_sample
    cat("Sample selected is: ", samp_name, "\n")
    samp_idx <<- which(samples() %in% samp_name)
    cat("Sample index is: ", samp_idx, "\n")
    output$data_viz = renderPlot({
      df_init = data()[data()$sample_ord==samp_name,]
      bar_plot_per_sample(data=df_init, sample_name = samp_name)
    })
    not_cols <<- c("genome_ord", "sample_name", "sample",
                   "run_id", "read", "sample_ord")
    if(is_single()==TRUE){
      cat("All reads are single-end. \n")
      output$r1_tbl = renderTable({
        df_init_r1 = num_data()[num_data()$sample_ord==samp_name,-which(colnames(num_data()) %in% not_cols)]
      })
      output$r2_tbl = renderTable({
        message = data.frame("Single-End reads. No R2 data to display.")
        colnames(message) = ""
        message
      })
    } else {
      output$r1_tbl = renderTable({
        df_init_r1 = num_data()[num_data()$sample_ord==samp_name & num_data()$read=="R1",]
        df_init_r1 = df_init_r1[,-which(colnames(df_init_r1) %in% not_cols)]
      })
      output$r2_tbl = renderTable({
        df_init_r2 = num_data()[num_data()$sample_ord==samp_name & num_data()$read=="R2",]
        df_init_r2 = df_init_r2[,-which(colnames(df_init_r2) %in% not_cols)]
      })
    }
  })
  
  # Move forward from current plot
  observeEvent(input$next_plt, {
    cat("------ NEXT PLOT BUTTON CLICKED ------", sep = "\n")
    if(samp_idx < samp_len){
      counter <<- samp_idx + 1
      cat("counter is set to: ", counter, "\n")
      samp_idx <<- counter
      cat("samp_idx is set to: ", samp_idx, "\n")
      output$data_viz = renderPlot({
        df_next = data()[data()$sample_ord==samples()[samp_idx],]
        bar_plot_per_sample(data=df_next, sample_name=samples()[samp_idx])
      })
      if(is_single()==TRUE){
        output$r1_tbl = renderTable({
          df_next_r1 = num_data()[num_data()$sample_ord==samples()[samp_idx],-which(colnames(num_data()) %in% not_cols)]
        })
        output$r2_tbl = renderTable({
          message = data.frame("Single-End reads. No R2 data to display.")
          colnames(message) = ""
          message
        })
      } else {
        output$r1_tbl = renderTable({
          df_next_r1 = num_data()[num_data()$sample_ord==samples()[samp_idx] & num_data()$read=="R1",]
          df_next_r1 = df_next_r1[,-which(colnames(df_next_r1) %in% not_cols)]
        })
        output$r2_tbl = renderTable({
          df_next_r2 = num_data()[num_data()$sample_ord==samples()[samp_idx] & num_data()$read=="R2",]
          df_next_r2 = df_next_r2[,-which(colnames(df_next_r2) %in% not_cols)]
        })
      }
    } else {
      cat("LAST PLOT VIEWED, RESETTING AT FIRST PLOT", sep = "\n")
      counter <<- 1
      cat("counter is set to: ", counter, "\n")
      samp_idx <<- counter
      cat("samp_idx is set to: ", samp_idx, "\n")
      output$data_viz = renderPlot({
        df_reset = data()[data()$sample_ord==samples()[samp_idx],]
        bar_plot_per_sample(data=df_reset, sample_name=samples()[samp_idx])
      })
      if(is_single()==TRUE){
        output$r1_tbl = renderTable({
          df_next_r1 = num_data()[num_data()$sample_ord==samples()[samp_idx],-which(colnames(num_data()) %in% not_cols)]
        })
        output$r2_tbl = renderTable({
          message = data.frame("Single-End reads. No R2 data to display.")
          colnames(message) = ""
          message
        })
      } else {
        output$r1_tbl = renderTable({
          df_next_r1 = num_data()[num_data()$sample_ord==samples()[samp_idx] & num_data()$read=="R1",]
          df_next_r1 = df_next_r1[,-which(colnames(df_next_r1) %in% not_cols)]
        })
        output$r2_tbl = renderTable({
          df_next_r2 = num_data()[num_data()$sample_ord==samples()[samp_idx] & num_data()$read=="R2",]
          df_next_r2 = df_next_r2[,-which(colnames(df_next_r2) %in% not_cols)]
        })
      }
    }
  })
  
  # Move backward from current plot
  observeEvent(input$last_plt, {
    cat("------ PREVIOUS PLOT BUTTON CLICKED ------", sep = "\n")
    if(samp_idx > 1){
      counter <<- samp_idx - 1
      cat("counter is set to: ", counter, "\n")
      samp_idx <<- counter
      cat("samp_idx is set to: ", samp_idx, "\n")
      output$data_viz = renderPlot({
        df_prev = data()[data()$sample_ord==samples()[samp_idx],]
        bar_plot_per_sample(data=df_prev, sample_name=samples()[samp_idx])
      })
      if(is_single()==TRUE){
        output$r1_tbl = renderTable({
          df_last_r1 = num_data()[num_data()$sample_ord==samples()[samp_idx],-which(colnames(num_data()) %in% not_cols)]
        })
        output$r2_tbl = renderTable({
          message = data.frame("Single-End reads. No R2 data to display.")
          colnames(message) = ""
          message
        })
      } else {
        output$r1_tbl = renderTable({
          df_last_r1 = num_data()[num_data()$sample_ord==samples()[samp_idx] & num_data()$read=="R1",]
          df_last_r1 = df_last_r1[,-which(colnames(df_last_r1) %in% not_cols)]
        })
        output$r2_tbl = renderTable({
          df_last_r2 = num_data()[num_data()$sample_ord==samples()[samp_idx] & num_data()$read=="R2",]
          df_last_r2 = df_last_r2[,-which(colnames(df_last_r2) %in% not_cols)]
        })
      }
    } else {
      cat("INITIAL PLOT VIEWED, MOVING BACKWARDS", sep = "\n")
      counter <<- samp_len
      cat("counter is set to: ", counter, "\n")
      samp_idx <<- counter
      cat("samp_idx is set to: ", samp_idx, "\n")
      output$data_viz = renderPlot({
        df_last = data()[data()$sample_ord==samples()[samp_idx],]
        bar_plot_per_sample(data=df_last, sample_name=samples()[samp_idx])
      })
      if(is_single()==TRUE){
        output$r1_tbl = renderTable({
          df_last_r1 = num_data()[num_data()$sample_ord==samples()[samp_idx],-which(colnames(num_data()) %in% not_cols)]
        })
        output$r2_tbl = renderTable({
          message = data.frame("Single-End reads. No R2 data to display.")
          colnames(message) = ""
          message
        })
      } else {
        output$r1_tbl = renderTable({
          df_last_r1 = num_data()[num_data()$sample_ord==samples()[samp_idx] & num_data()$read=="R1",]
          df_last_r1 = df_last_r1[,-which(colnames(df_last_r1) %in% not_cols)]
        })
        output$r2_tbl = renderTable({
          df_last_r2 = num_data()[num_data()$sample_ord==samples()[samp_idx] & num_data()$read=="R2",]
          df_last_r2 = df_last_r2[,-which(colnames(df_last_r2) %in% not_cols)]
        })
      }
    }
  })

  # Render UI outputs
  output$tb = renderUI({
    tabsetPanel(tabPanel("Uploaded Files", tableOutput("files_df")),
                tabPanel("Percent Mapping Data",
                         DT::dataTableOutput("data_summary"),
                         br(),
                         downloadButton('data_dl',"Download this table"),
                         tags$style(type='text/css', "#data_dl { color: #00A1E6; }")),
                tabPanel("Plots",
                         selectInput("which_sample", label = "Select a sample:",
                                     choices = data()$sample_ord),
                         plotOutput("data_viz", height = "700px"),
                         actionButton("next_plt", "Next Sample >>"),
                         actionButton("last_plt", "<< Previous Sample"),
                         tabsetPanel(tabPanel("R1 Number Mapping",
                                              tableOutput("r1_tbl")),
                                     tabPanel("R2 Number Mapping",
                                              tableOutput("r2_tbl")))))
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)