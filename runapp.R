shiny::runApp(appDir = getwd(), 
       launch.browser = getOption("shiny.launch.browser", interactive()),
       host = getOption("shiny.host", "172.21.162.112"), workerId = "",
       quiet = FALSE, display.mode = c("auto", "normal", "showcase"),
       test.mode = getOption("shiny.testmode", FALSE))
