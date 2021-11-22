shiny::runApp(appDir = getwd(), port = 9000,
       launch.browser = getOption("shiny.launch.browser", interactive()),
       host = getOption("shiny.host", "0.0.0.0"), workerId = "",
       quiet = FALSE, display.mode = c("auto", "normal", "showcase"),
       test.mode = getOption("shiny.testmode", FALSE))


