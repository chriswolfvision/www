#!/usr/bin/ruby
# ********************************************************************
# lib.rb
# A ruby library of handy helper functions
#
# Author: Christian Wolf, christian.wolf@insa-lyon.fr
# ********************************************************************

# **************************************************************************
# A couple of definitions
# **************************************************************************

$RubyPath="/usr/bin/ruby"

# **************************************************************************
# Run the command given in the string
# **************************************************************************

def runCom(com,doEcho)
	if doEcho 
		$stderr.printf "%s\n",com 
		$stderr.flush
	end
	if not system(com)		
		$stderr.printf "Error code: %d\n",$?				
		exit 1
	end
end

# **************************************************************************
# returns a filename w/o it's suffix
# **************************************************************************

def noSuffix(s)
	return s[0..s.rindex(".")-1]
end

# **************************************************************************
# Array 2 string
# **************************************************************************

def arr2str(arr)
	rv=""
	arr.each do |i|
		rv=rv+i+" "
	end
	return rv
end


