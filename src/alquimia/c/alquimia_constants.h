/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_STRINGS_H_
#define ALQUIMIA_STRINGS_H_

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /* String lengths */
  static const int ALQUIMIA_MAX_STRING_LENGTH = 512;
  static const int ALQUIMIA_MAX_WORD_LENGTH = 32;

  /* Geochemistry Engine Strings */
  static const char* kPFloTran = "PFloTran";
  static const char* kCrunchFlow = "CrunchFlow";

  /* Error Codes */
  static const int ALQUIMIA_NO_ERROR = 0;
  static const int ALQUIMIA_ERROR_INVALID_ENGINE = 1;
  static const int ALQUIMIA_ERROR_ENGINE_INTEGRITY = 4577;

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif     /* ALQUIMIA_STRINGS_H_ */
