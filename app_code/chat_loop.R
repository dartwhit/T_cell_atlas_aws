run_tool_loop <- function(messages, system_prompt, seurat_obj, tool_schemas,
                           max_iterations = 8,
                           model = "claude-haiku-4-5-20251001") {
  api_key <- Sys.getenv("ANTHROPIC_API_KEY")
  if (nchar(api_key) == 0) stop("ANTHROPIC_API_KEY environment variable is not set.")

  MAX_HISTORY_TURNS <- 20
  current_messages <- if (length(messages) > MAX_HISTORY_TURNS) {
    tail(messages, MAX_HISTORY_TURNS)
  } else {
    messages
  }

  for (i in seq_len(max_iterations)) {
    body <- list(
      model      = model,
      max_tokens = 2048,
      system     = system_prompt,
      tools      = tool_schemas,
      messages   = current_messages
    )

    resp <- tryCatch(
      httr2::request("https://api.anthropic.com/v1/messages") |>
        httr2::req_headers(
          "x-api-key"         = api_key,
          "anthropic-version" = "2023-06-01",
          "content-type"      = "application/json"
        ) |>
        httr2::req_body_json(body) |>
        httr2::req_timeout(120) |>
        httr2::req_perform(),
      error = function(e) {
        stop(paste("API request failed:", conditionMessage(e)))
      }
    )

    result <- httr2::resp_body_json(resp)

    if (!is.null(result$error)) {
      return(paste0("**API error:** ", result$error$message))
    }

    if (result$stop_reason == "end_turn") {
      text_blocks <- Filter(function(b) b$type == "text", result$content)
      return(paste(sapply(text_blocks, `[[`, "text"), collapse = "\n"))
    }

    if (result$stop_reason == "tool_use") {
      current_messages <- c(current_messages, list(list(
        role    = "assistant",
        content = result$content
      )))

      tool_results <- lapply(
        Filter(function(b) b$type == "tool_use", result$content),
        function(block) {
          output <- dispatch_tool(block$name, block$input, seurat_obj)
          list(
            type        = "tool_result",
            tool_use_id = block$id,
            content     = output
          )
        }
      )

      current_messages <- c(current_messages, list(list(
        role    = "user",
        content = tool_results
      )))
      next
    }

    break
  }

  "[Error: tool loop exceeded maximum iterations or received an unexpected stop reason]"
}
