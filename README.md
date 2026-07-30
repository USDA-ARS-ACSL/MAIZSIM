# MAIZSIM

A mechanistic model of maize growth, development, and yield  
Developed by the USDA‑ARS Adaptive Cropping Systems Laboratory and the University of Washington School of Environmental and Forest Sciences

***

## Overview

MAIZSIM is a mechanistic model of maize growth, development, and yield.  
The crop system is written in **C++**, and the soil system is written in **Fortran**.  
The crop model is integrated with **2DSOIL**, a two‑dimensional simulator of soil water and heat movement and solute transport.  
**2DSOIL is the main driver model**, calling the crop model as a subroutine.

The repository includes two subprojects:

* **Crop Source**
* **Soil Source**

***

## Developers

### Primary Developers

* Soo‑Hyung Kim — University of Washington
* Dennis Timlin — USDA‑ARS
* David Fleisher — USDA‑ARS
* V. R. Reddy — USDA‑ARS
* Zhuangji Wang - Colorado State Univesity

### Additional Contributors

* Yang Yang — Dow Agrosciences
* Annette Dathe — Cornell University *(diffusive root model)*
* Jong‑Ahn Chun — APEC Climate Center, Korea *(CO₂ and water)*
* Sahila Beegum — University of Nebraska *(gas transport & respiration)*
* Wenguang Sun — Colorado State University *(gas transport & respiration)*

The MAIZSIM community encourages **responsible disclosure of vulnerabilities** to improve security and reliability.

***

## Branch Structure

* **dev** — holds recent code updates
* **master** — stable version  
  Changes are tested in **dev** and merged into **master** after verification.

**To collaborate:**

1. Fork the repository
2. Create a feature branch
3. Submit a pull request

***

## Building the Model

### Windows Build (Visual Studio + Intel Fortran)

MAIZSIM builds successfully using:

* **Visual Studio Professional 2022**
* **Intel Fortran OneAPI‑2023.2**

Macros copy the compiled crop DLLs into the soil project's directory.  
To ensure proper library linking, the directory structure must remain:

```
crop source/
soil source/
```

### Linux Build (Docker)

Linux compatibility was added by **Kyungdahm Yun** (University of Washington).  
A Docker image containing compilers and makefiles is available:

👉 <https://github.com/precision-sustainable-ag/BuildMaizsim>

This image may be used under Linux or through **WSL2** on Windows.

***

## Running the Model

More documentation is in preparation.  
See **how to run the model** in the repository for details on:

* Preparing input files
* Running the executable
* Command‑line usage

***

## Excel Interface

An Excel‑based interface, including example input files, is available here:

👉 <https://github.com/USDA-ARS-ACSL/ExcelInterface>

This interface supports the most recent MAIZSIM version.

***

## Updates Since 2025

| Date                      | Description                                                                               |
| ------------------------- | ----------------------------------------------------------------------------------------- |
| **7/7/2026 12:17:33 PM**  | Increased max fertilizer applications to 175                                              |
| **3/23/2026 3:11:46 PM**  | Updated version number in splash screens                                                  |
| **3/20/2026 9:52:28 AM**  | Redid NH₄ changes (previous updates pulled from incorrect repo)                           |
| **3/20/2026 9:23:31 AM**  | Fixed NH₄ initialization error caused by BD not yet read                                  |
| **9/24/2025 11:59:52 AM** | Added comment                                                                             |
| **9/24/2025 11:59:14 AM** | Harmonized single/double precision; fixed NaN from `dt` being double; cleaned unused vars |
| **8/12/2025 5:09:25 PM**  | Added safeguard for excessive daily infiltration                                          |
| **7/23/2025 4:10:22 PM**  | Linux compatibility fixes; removed outdated Fortran constructs                            |
| **5/13/2025 9:53:02 AM**  | Added date to header                                                                      |
| **5/13/2025 9:52:26 AM**  | Added double descriptor to `abio` call for Linux                                          |
| **5/13/2025 9:51:05 AM**  | Renamed weather structure to avoid Linux conflicts                                        |
| **5/13/2025 9:33:46 AM**  | Updated weather common block name                                                         |
| **3/7/2025 9:13:15 AM**   | Added more comments                                                                       |
| **1/15/2025 1:18:14 PM**  | Added constraints for Tmin/Tmax when below 0°C                                            |

***

## License

*CC0 1.0 Universal (Public Domain Dedication)
--------------------------------------------

To the extent possible under law, the authors have waived all copyright
and related or neighboring rights to MAIZSIM. This work is published
from the United States.

You may copy, modify, distribute, and perform the work, even for
commercial purposes, all without asking permission.
For full details, see: https://creativecommons.org/publicdomain/zero/1.0/*

***

## Vulnerability Disclosure Policy
**The community is explicitly encouraged to engage in the responsible disclosure of vulnerabilities to promote collaboration and improve code security.**
If you discover a security vulnerability, please report it responsibly:
1.	**Email:** name.lastname@usda.gov with subject "Security Vulnerability - Probit Tool"
1.	**Provide details:** Description, steps to reproduce, potential impact
1.	**Confidential handling:** We will respond within 48 hours
1.	**Recognition:** Contributors acknowledged (with permission) after resolution

**Please do not publicly disclose vulnerabilities until they have been addressed.**
**Vulnerability Response Timeline**
When vulnerabilities are identified:
+	**Critical vulnerabilities:** Patched within 7 days or application taken offline
+	**High vulnerabilities:** Addressed within 14 days
+	**Medium/Low vulnerabilities:** Resolved within 30 days
+	**Users notified:** Via GitHub releases and repository notices
+	**Workarounds provided:** If immediate fixes are not possible
**If vulnerabilities cannot be timely resolved, a prominent warning will be added to this README and the application may be temporarily taken offline until fixes are implemented.**

## Contact

For questions, improvements, and collaboration, please open an issue or pull request on GitHub.

***


