library(shiny)

source("ui.R")
source("server.R")

# Wrap UI with shinymanager authentication (unless disabled via environment variable)
disable_auth <- Sys.getenv("DISABLE_AUTH", "false") == "true"
if (disable_auth) {
  cat("⚠️  Authentication DISABLED via DISABLE_AUTH environment variable\n", file = stderr())
  # ui is already defined from ui.R, no wrapping needed
} else {
  cat("✅ Authentication ENABLED - wrapping UI with secure_app()\n", file = stderr())
  library(shinymanager)
  
  # Extract custom head tags before wrapping
  custom_scripts <- tags$head(
    tags$script(HTML('
      // Wait for jQuery to be available before using it
      function waitForJQuery(callback) {
        if (typeof jQuery !== "undefined") {
          callback(jQuery);
        } else {
          setTimeout(function() { waitForJQuery(callback); }, 50);
        }
      }
      
      // Use jQuery once it\'s available
      waitForJQuery(function($) {
        $(document).on("shiny:connected", function() {
          function updateDimensions() {
            Shiny.onInputChange("dimension", [window.innerWidth, window.innerHeight]);
          }
          updateDimensions();
          $(window).resize(updateDimensions);

          // Any element with id starting with "explore_" switches to Explore tab
          $("[id^=\'explore_\']").on("click", function() {
            $("a[data-value=\'Explore\']").tab("show");
          });
        });
      });
    '))
  )
  
  tryCatch({
    # Use tags_top parameter to ensure custom scripts are included
    ui <- secure_app(ui, enable_admin = TRUE, tags_top = custom_scripts)
    cat("✓ UI wrapped with authentication successfully\n", file = stderr())
  }, error = function(e) {
    cat("❌ ERROR wrapping UI with secure_app():", conditionMessage(e), "\n", file = stderr())
    cat("⚠️  Continuing with unwrapped UI\n", file = stderr())
    # ui remains as the original page_navbar() object
  })
}

# Verify ui is defined
if (!exists("ui")) {
  stop("ERROR: ui object was not found after sourcing ui.R!")
}
cat("✓ UI object ready for shinyApp\n", file = stderr())

shinyApp(ui = ui, server = server)