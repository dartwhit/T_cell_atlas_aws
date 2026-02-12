# Disabling Authentication for PR/Testing Deployments

For temporary PR deployments or local testing, you can disable the shinymanager login authentication using the `DISABLE_AUTH` environment variable.

## Quick Start

### Option 1: Using docker-compose.override.yml (Recommended)

1. Copy the example override file:
   ```bash
   cp docker-compose.override.yml.example docker-compose.override.yml
   ```

2. Run docker-compose as normal:
   ```bash
   docker-compose up -d
   ```

The override file will automatically disable authentication for all shiny app instances.

### Option 2: Setting environment variable directly

Add the environment variable to your docker-compose.yml or set it when running:

```yaml
services:
  shiny-app-1:
    environment:
      - DISABLE_AUTH=true
```

### Option 3: For local R development

Set the environment variable before running the app:

```r
Sys.setenv(DISABLE_AUTH = "true")
shiny::runApp()
```

Or in your R session:
```bash
export DISABLE_AUTH=true
R -e "shiny::runApp('app_code')"
```

## Verification

When authentication is disabled, you should see in the logs:
```
⚠️ Authentication disabled via DISABLE_AUTH environment variable
```

And the app will load directly without showing the login screen.

## Security Note

⚠️ **Never disable authentication in production deployments!** This feature is intended only for:
- PR preview deployments
- Local development
- Testing environments

For production, always ensure `DISABLE_AUTH` is not set or is set to `false`.
