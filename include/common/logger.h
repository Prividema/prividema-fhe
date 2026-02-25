#ifndef LOGGER_H
#define LOGGER_H

/**
 *  @brief Enum that defines the different log levels :
 * `LOG_DEBUG`, `LOG_INFO`, `LOG_WARN`, `LOG_ERROR`
 */
typedef enum { LOG_DEBUG, LOG_INFO, LOG_WARN, LOG_ERROR } LogLevel;

extern LogLevel LOG_THRESHOLD;

/**
 * @brief Set the program's log level.
 * @param level The log level. Should be an element of the enum.
 */
void set_log_level(LogLevel level);

/**
 * @brief Log a message with string formatting.
 * @param level The log level.
 * @param fmt A va_list with the string to format and its arguments.
 * @retval -1 if the level is LOG_ERROR.
 * @retval  0 otherwise.
 */
int log_msg(LogLevel level, const char* fmt, ...);

/**
 * @brief Log a syscall error like perror but with LOG_ERROR.
 * @param header The perror message header.
 * @return always -1 as a syscall error is often fatal.
 */
int log_perror(const char* header);

/**
 * @brief Log a syscall error like log_perror, if and only if ptr is NULL.
 * @param ptr The pointer to test.
 * @param header The perror message header.
 * @return -1 if and only if test is NULL.
 */
int log_is_null(void* ptr, const char* header);
#endif  // LOGGER_H