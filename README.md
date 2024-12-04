# Seq2VCF Pipeline

Le **Seq2VCF Pipeline** est un outil bioinformatique conçu pour traiter des données de séquençage brut et générer des fichiers Variant Call Format (VCF) de haute qualité. Ce pipeline automatise les étapes nécessaires, telles que la démultiplexation, le contrôle qualité, l’alignement, l’appel de variants et l’annotation fonctionnelle.

---

## Objectif et Fonctionnalités

Le pipeline a été développé pour traiter efficacement les données de génotypage par séquençage (GBS), permettant l’identification rapide et reproductible des variants génétiques. Ses principales caractéristiques sont :  
- Traitement automatisé des lectures brutes jusqu’aux variants annotés.  
- Prise en charge de plusieurs échantillons via la démultiplexation des codes-barres.  
- Contrôle qualité complet et gestion des erreurs.

---

## Formats des Fichiers d'Entrée et de Sortie

### Formats d’entrée :
1. **Données de séquençage brutes** : Fichiers FASTQ compressés (`.fq.gz`) contenant les lectures.  
2. **Fichier de codes-barres** : Fichier texte tabulé contenant les codes-barres des échantillons (`barcodes.txt`).  
3. **Génome de référence** : Fichier FASTA (`.fa`) contenant la séquence du génome.

### Formats de sortie :
1. **Lectures démultiplexées** : Fichiers FASTQ (`.fq.gz`) par échantillon.  
2. **Rapports qualité** : Rapports HTML et texte générés par FastQC.  
3. **Lectures alignées** : Fichiers BAM (`.bam`) par échantillon.  
4. **Variants détectés** : Fichiers VCF (`.vcf`) contenant les variants détectés et filtrés.  
5. **Variants annotés** : Fichiers VCF annotés (`.annotated.vcf`) avec des informations fonctionnelles (si activé).

---

## Dépendances et Étapes d’Installation

### Logiciels requis :
Le pipeline nécessite les outils suivants :  
- **Python** (v3.5+)  
- **Java** (v22.0.2+)  
- **Cutadapt** (v2.1+)  
- **Sabre** (v1.000+)  
- **BWA** (v0.7.17+)  
- **SAMtools** (v1.8+)  
- **BCFtools** (v1.15+)  
- **FastQC** (v0.11.2+)  
- **SnpEff** (v5.2e)

### Étapes d’installation :
1. Cloner le dépôt :
   ```bash
   git clone https://github.com/username/Seq2VCF.git
   cd Seq2VCF
   ```
2. if necessary, install required tools.  Use conda, apt, or brew based on your system. Example for conda:
     ```bash
   conda install -c bioconda bwa samtools bcftools fastqc cutadapt
   ```
3. Ensure the tools are accessible in your PATH by adding them to your .bashrc or .zshrc.

---
## How to Run the Pipeline
### Configuration:
Edit the **parameters.conf** file to specify input paths, output directories, and processing options. Key parameters include:

**DATA_DIR**: Path to the raw data folder.
**BARCODES_FILE**: Path to the barcode file.
**REF_GENOME**: Path to the reference genome file.
**RESULTS_DIR**: Directory for saving results.
**RUN_ANNOTATION**: Enable functional annotation (true/false).
### Execution:
Run the pipeline using the main script:

```bash
Copier le code
bash pipeline.sh
```
---
## Log File Details and Troubleshooting
The pipeline generates a log file (pipeline.log) for tracking progress and debugging.

### Log File Contents:
**INFO**: Normal operation messages.
**WARN**: Non-critical issues that may need attention.
**ERROR**: Critical errors that terminate the pipeline.
### Troubleshooting:
**Missing Dependencies**: Ensure all tools are installed and accessible in the PATH.
**File Not Found**: Check the parameters.conf file for correct paths.
**Annotation Errors**: If SnpEff fails, ensure the correct genome database is downloaded and specified in parameters.conf.
**Out of Memory**: Allocate sufficient memory to Java-based tools (e.g., SnpEff). Edit JAVA_OPTS in parameters.conf if needed.

For additional help, refer to the documentation in the **docs/** folder or contact the maintainer.
