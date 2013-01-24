These are the named constants that alquimia uses and their required values. 

== **Alquimia Error Codes** ==
|=name |=value |=meaning
|ALQUIMIA_NO_ERROR | 0 | no error
|ALQUIMIA_ERROR_ENGINE_INTEGRITY | 4577 | pointer to the engine's internal state did not pass integrity check

== **Alquimia Strings** ==

|=name |=value
| kPFloTran | "PFloTran"
| kCrunchFlow | "CrunchFlow"
| kCppChem | "C++ Chem"
| kpH | "pH"
| kMineral | "mineral"
| kGas | "gas"
| kCharge | "charge"

NOTE(bja): Except for species names, strings should be case insensitive? The demo driver (and C++ chem) use lower case, but pflotran uses a mixture of StringToUpper and StringToLower internally.

== **Alquimia String Lengths**==
  static const int ALQUIMIA_MAX_STRING_LENGTH = 512;
  static const int ALQUIMIA_MAX_WORD_LENGTH = 32;
