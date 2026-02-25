#include "logger.h"

#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

// Global log level threshold
LogLevel LOG_THRESHOLD = LOG_DEBUG;

// log threshold setter
void set_log_level(LogLevel level) { LOG_THRESHOLD = level; }

// Get string representation of log level
const char* log_level_str(LogLevel level)
{
	switch (level) {
		case LOG_DEBUG:
			return "DEBUG";
		case LOG_INFO:
			return "INFO";
		case LOG_WARN:
			return "WARN";
		case LOG_ERROR:
			return "ERROR";
		default:
			return "UNKNOWN";
	}
}

// The main logging function
int log_msg(LogLevel level, const char* fmt, ...)
{
	// Skip lower-priority messages
	if (level < LOG_THRESHOLD) return 0;

	// Get current time
	time_t t = time(NULL);
	struct tm tm_info;
	localtime_r(&t, &tm_info);  // thread-safe version

	char time_buf[20];
	strftime(time_buf, sizeof(time_buf), "%Y-%m-%d %H:%M:%S", &tm_info);

	// Print log level and timestamp
	printf("[%s] %s: ", time_buf, log_level_str(level));

	// Handle the variable arguments
	va_list args;
	va_start(args, fmt);
	vprintf(fmt, args);
	va_end(args);

	printf("\n");  // newline at the end

	if (level == LOG_ERROR) return -1;
	return 0;
}

// Perror with logger
int log_perror(const char* header) { return log_msg(LOG_ERROR, "%s : %s", header, strerror(errno)); }

int log_is_null(void* ptr, const char* header)
{
	if (ptr == NULL) {
		return log_perror(header);
	} else
		return 0;
}