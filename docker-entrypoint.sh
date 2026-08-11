#!/bin/sh
# Shiny Server spawns the R worker with a hardcoded environment whitelist
# (HOME, LANG, PATH only), so anything passed via `docker run -e` or the
# compose `environment:` block is invisible to Sys.getenv() inside the app.
# Materialize the vars we care about into Renviron.site, which R reads on
# every startup, then hand off to shiny-server.
set -eu

APP_ENV_VARS="LOCAL_DEV ANTHROPIC_API_KEY QS2_READ_THREADS"

RENVIRON_SITE="$(R RHOME)/etc/Renviron.site"
touch "$RENVIRON_SITE"

for var in $APP_ENV_VARS; do
    eval "value=\${$var:-}"
    [ -n "$value" ] || continue
    # Drop any previous definition so restarts don't accumulate duplicates.
    grep -v "^${var}=" "$RENVIRON_SITE" > "${RENVIRON_SITE}.tmp" || true
    mv "${RENVIRON_SITE}.tmp" "$RENVIRON_SITE"
    printf '%s=%s\n' "$var" "$value" >> "$RENVIRON_SITE"
    echo "entrypoint: exported ${var} to ${RENVIRON_SITE}"
done

exec "$@"
