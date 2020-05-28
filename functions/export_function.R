csv_export_handler<-function(data){downloadHandler(
  filename = function() {
    'pat.csv'
  },
  content = function(file) {
    write.csv(data, file)
  }
)
}


