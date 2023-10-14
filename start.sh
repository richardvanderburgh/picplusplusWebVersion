#!/bin/sh

# Check if $PORT is set, if not, set to 8000
if [ -z "$PORT" ]; then
    export PORT=8000
fi

echo "Port is: $PORT"

# Start Gunicorn
exec gunicorn djangodocker.wsgi:application --bind 0.0.0.0:$PORT
