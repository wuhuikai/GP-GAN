#!/bin/bash
cat >._lighttpd_dir_conf << EOF
server.document-root = "$1"

server.port = "$2"

dir-listing.encoding        = "utf-8"
server.dir-listing          = "enable"

mimetype.assign = (
  ".html" => "text/html",
  ".txt" => "text/plain",
  ".jpg" => "image/jpeg",
  ".png" => "image/png"
)
EOF

lighttpd -D -f ._lighttpd_dir_conf