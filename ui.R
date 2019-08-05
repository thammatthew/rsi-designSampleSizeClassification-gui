dashboardPage(
  skin="black",
  dashboardHeader(
    title="designSampleSizeClassification",
    titleWidth=300,
    tags$li(class = "dropdown",
            tags$li(class = "dropdown", actionLink("debug_save", label="Save Environment"))
    )
  ),

  dashboardSidebar(
    width=300,
    sidebarMenu(id="tabs",
      menuItem("Home", tabName="home", icon=icon("home")),
      menuItem("Import Data", icon=icon("file-import"), startExpanded=TRUE,
        div(class="custom-sidebar-output",
          selectInput("data_format", "Select Data Format", choices=list("Annotated Protein-level Data"="standard", "MSstats Data"="MSstats", "Examples (MSstatsBioData)"="Examples", "Debug"="Debug"))
        ),
        menuItem("Options", icon=icon("sliders-h"),
          htmlOutput("options", class="custom-sidebar-output")
        ),
        menuItem("Advanced", icon=icon("cogs"),
          htmlOutput("advanced_options", class="custom-sidebar-output")
        ),
        htmlOutput("select_files", class="custom-sidebar-output",)
      ),
      menuItem("Explore Raw Data", tabName="explore_data", icon = icon("search")),
      menuItem("Explore Simulated Data", tabName="explore_simulated", icon=icon("project-diagram")),
      menuItem("Analyse Simulated Data", tabName="analyse_simulated", icon=icon("vial"))
    )
  ),
  
  dashboardBody(
    useShinyjs(),
    tags$head(
      tags$style(HTML("
      .box {
        overflow-y: scroll;
      }
      /* Styling for elements nested inside menuItems */
      .treeview-menu .shiny-input-container {
        margin-bottom:0px !important;
        padding-top:0px !important;
        padding-bottom:0px !important;
      }
      .custom-sidebar-output>.shiny-input-container .progress{
        margin-bottom:0px !important;
      }
      .custom-sidebar-output>b {
        padding-left:15px;
      }
      .custom-sidebar-output {
        margin: 10px 0 10px 0;
      }
      .treeview-menu {
        padding:5px 0px 5px 5px !important;
      }

      #main_output {
        font-family: Monaco, Consolas, 'Andale Mono', 'DejaVu Sans Mono', monospace;
        font-size: 80%;
        margin-top: 5px;
        white-space: pre;
        white-space: pre-wrap;
        white-space: -moz-pre-wrap;
        white-space: -o-pre-wrap;
        background: #222D32;
        color: #FFFFFF;
        border-radius: 2px;
      }
      
      .content-wrapper, .right-side {
        background-color: #F7F7F7;
      }
      
      h1 {
        margin-top: 0px;
      }
      
      * {
        word-break: break-word;
      }
      .custom-box-title {
        display: inline-block;
        font-size: 18px;
        margin: 0;
        line-height: 1;
      }
    "))),
    tabItems(
      Tab1 <- tabItem(
        tabName="home",
        includeMarkdown("www/Welcome.md")
      ),
      
      Tab2 <- tabItem(
        tabName="explore_data",
        htmlOutput("explore_data_content")
      ),
      
      Tab3 <- tabItem(
        tabName="explore_simulated",
        h1("Explore Simulated Datasets"),
        htmlOutput("explore_simulated_content")
      ),
      
      Tab4 <- tabItem(
        tabName = "analyse_simulated",
        h1("Analyse Simulated Data"),
        htmlOutput("analyse_simulated_content")
      )
    )
  )
)