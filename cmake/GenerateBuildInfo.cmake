# Queries git for the current commit and working-tree state, then expands
# build_info.hpp.in. Run at build time (not configure time) so the recorded
# commit does not go stale after every commit; configure_file only rewrites
# the output when its content actually changes, so this does not force a
# rebuild on every invocation.
#
# Degrades to "unknown" when git is missing or the tree is not a repository:
# a missing commit id is not a reason to fail someone's build.

find_package(Git QUIET)

set(GIT_COMMIT "unknown")
set(GIT_DIRTY "")

if(GIT_FOUND)
	execute_process(
		COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
		WORKING_DIRECTORY ${GIT_DIR}
		OUTPUT_VARIABLE GIT_REV
		OUTPUT_STRIP_TRAILING_WHITESPACE
		RESULT_VARIABLE GIT_REV_RESULT
		ERROR_QUIET
	)
	if(GIT_REV_RESULT EQUAL 0 AND GIT_REV)
		set(GIT_COMMIT "${GIT_REV}")

		execute_process(
			COMMAND ${GIT_EXECUTABLE} status --porcelain
			WORKING_DIRECTORY ${GIT_DIR}
			OUTPUT_VARIABLE GIT_STATUS
			OUTPUT_STRIP_TRAILING_WHITESPACE
			ERROR_QUIET
		)
		if(GIT_STATUS)
			set(GIT_DIRTY "-dirty")
		endif()
	endif()
endif()

configure_file(${SRC} ${DST} @ONLY)
