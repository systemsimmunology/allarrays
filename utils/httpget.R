

hex.value <- function( s ) {
	as.integer( paste("0x", s, sep="" ) );
}

# This handles 2 of 4 possible entity-body cases:
# 1) Content-Length header indicates size of entity-body
# 2) Transfer-Encoding: chunked indicates a chunked entity-body
# Note that this code pays attention to neither chunked entity body
# trailers NOR the transfer-extension.

# Pseudo-code for receiving chunked data:
# (from http://www.w3.org/Protocols/rfc2616/rfc2616-sec19.html#sec19.4.6) 
#
# length := 0
# read chunk-size, chunk-extension (if any) and CRLF
# while (chunk-size > 0) {
#   read chunk-data and CRLF
#   append chunk-data to entity-body
#   length := length + chunk-size
#   read chunk-size and CRLF
# }
# read entity-header
# while (entity-header not empty) {
#   append entity-header to existing header fields
#   read entity-header
# }
# Content-Length := length
# Remove "chunked" from Transfer-Encoding

http.get <- function( host, path, port=80, debug=FALSE ) {

	# Might want to provide a timeout optional argument, and do...
	#   .Options$timeout <- timeout;
	# ...here.

	if( missing(path) ) path <- "/";

	fp  <- socketConnection( host=host, port=port, server=FALSE, blocking=TRUE );
	stopifnot( isOpen( fp ) );

	# Build and write the request header

	hdr <- character(0);
	hdr <- c(hdr, paste("GET ",path," HTTP/1.1\n",sep=""))
	hdr <- c(hdr, paste("Host: ",host,"\n",sep=""))
	hdr <- c(hdr, "User-Agent: R http.get function\n")
	hdr <- c(hdr, "Accept: */*\n\n")
	hdr <- paste( hdr, collapse="" )

	if( debug ) {
		cat( "+++ Request Headers...\n" );
		cat( "<<< ", hdr );
	}
	writeChar( hdr, fp, nchars=nchar(hdr,type="chars"), eos=NULL );

	# Read the response...
	# ...first the headers which, happily, are line-based.

	entity.header <- character(0);

	repeat {
		l <- readLines( con=fp, 1, ok=FALSE );
		if( debug ) cat( ">>> ", l, "\n" );
		# NOTE! readLines appears to consume CRNL pairs like the C function gets().
		entity.header <- c( entity.header, l );
		if( nchar(l) == 0 || '\r\n' == l ) {
			break; # 1st empty line terminates header
		}
	}

	content.length <- integer(1);

	for( l in entity.header ) {
		# Look for the Content-Length or Transfer-Encoding header lines.
		m <- regexpr( "Content-Length: *[0-9]+", l );
		if( m > -1 ) {
			# ...extract the content length value.
			m.len <- attr( m, 'match.length' );
			i <- m.len;
			while( regexpr( '[0-9]', substr(l,i,i) ) > -1 ) { 
				i <- i-1 
			}
			content.length[[1]] <- as.integer( substr(l, i+1, m.len ) );
			next;
		}
		m <- regexpr( "Transfer-Encoding: *chunked", l );
		if( m > -1 ) {
			content.length[[1]] <- -1;
			# "If a message is received with both a Transfer-Encoding header 
			# field and a Content-Length header field, the latter MUST be 
			# ignored," so...
			break;
		}
	}

	if( debug ) cat( ">>> beginning read of entity.body...\n" );

	entity.body <- if( content.length[[1]] < 0 ) {

		if( debug ) cat( ">>> ...in chunked transfer\n" );

		# Read chunked data.
		chunked.body <- character(0);

		repeat {
			l <- readLines( con=fp, 1, ok=FALSE ); # chunk-size [ chunk-extension ] and CRLF
			m <- regexpr( '[0-9a-fA-F]+', l );
			clen <- hex.value( substr( l, m, attr( m, 'match.length' ) ) );
			if( clen == 0 ) 
				break;
			chunked.body <- c( chunked.body, readChar( fp, clen ) );
			l <- readLines( con=fp, 1, ok=FALSE ); # CRLF
		}

		# Read the optional trailer.
		# Actually I'm just consuming the trailer.

		if( debug ) cat( ">>> beginning read of entity.body trailer\n" );

		repeat {
			l <- readLines( con=fp, 1, ok=FALSE );
			if( nchar(l) > 0 ) {
				# Careful with braces and semicolons here! They affect
				# control flow in surprising ways!
				if( debug ) cat( ">>> ", l, "\n" );
			} else
				break;
		}

		paste( chunked.body, collapse="" );

	} else {

		if( debug ) cat( ">>> ...with Content-Length", content.length[[1]], "\n" );
		readChar( fp, content.length[[1]] );
	}

	close( fp );
	return(entity.body);
}

