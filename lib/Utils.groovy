import groovy.json.JsonSlurper
import groovy.csv.CsvParser

class Utils {
    static def jsonSlurper = new JsonSlurper()

    static def checkResources(config) {
        def globresource = "bioit"
        if (config.globalResources.contains("references_backup") && !config.containsKey("organism")) {
            if (config.containsKey("reference") && config.reference.contains("_r")) {
                globresource = "bioit"
                config.globalResources = config.globalResources.replace("base/references_backup", "resources")
            } else {
                globresource = "bioda"
            }
        } else {
            globresource = "bioit"
            config.globalResources = config.globalResources.replace("base/references_backup", "resources")
        }
        return config + [globresource: globresource]
    }

    static def loadSample(config) {
        return config.samples.collectEntries { k, v -> [k, v] }
    }

    static def setReadPairTags(config) {
        return config.is_paired ? ["_R1", "_R2"] : [""]
    }

    static def setReadPairQcTags(config) {
        return config.is_paired ? ["R1", "R2"] : ["SE"]
    }

    static def setPairedTags(config) {
        return config.is_paired ? "PE" : "SE"
    }

    static def setReadPairDmtexTags(config) {
        return config.is_paired ? ["_R1", "_R2"] : ["_R1"]
    }

    static def loadRef(config) {
        if (config.lib_ROI != "wgs") {
            def libRoiDict = jsonSlurper.parseText(new File("${config.globalResources}/reference_info/lib_ROI.json").text)
            config.reference = libRoiDict.find { k, v -> v instanceof Map && v.containsKey(config.lib_ROI) }?.key
        }
        return config
    }

    static def loadTooldir(config) {
        def globresource = config.globresource
        if (globresource == "bioda") {
            config.tooldir = "${config.globalResources}/general"
        } else if (globresource == "bioit") {
            config.tooldir = "${config.globalResources}/tools"
        }
        return config
    }

    static def loadROI(config) {
        def libRoiDict = jsonSlurper.parseText(new File("${config.globalResources}/reference_info/lib_ROI.json").text)
        
        if (config.lib_ROI == "rna") {
            config.material = "RNA"
        } else {
            if (config.globresource == "bioda") {
                config.reference = libRoiDict.find { k, v -> v instanceof Map && v.containsKey(config.lib_ROI) }?.key
            }
            if (config.globresource == "bioit") {
                if (!config.containsKey("organism")) {
                    config.reference = libRoiDict.find { k, v -> v instanceof Map && v.containsKey(config.lib_ROI) }?.key
                }
                config.lib_ROI = config.lib_ROI.split("_")[0..-2].join("_")
            }
        }
        return config
    }

    static def loadOrganism(config) {
        def referenceDict = jsonSlurper.parseText(new File("${config.globalResources}/reference_info/reference2.json").text)
        def keggDict = jsonSlurper.parseText(new File("${config.globalResources}/reference_info/kegg_reference.json").text)
        def organismTab = new File("${config.globalResources}/reference_info/organism_tab.tsv").readLines()
            .drop(1) // skip header
            .collect { line -> 
                def (assembly, full_name, kegg_term, release) = line.split('\t')
                [assembly: assembly, full_name: full_name, kegg_term: kegg_term, release: release]
            }
            .groupBy { it.assembly }
            .collectEntries { assembly, list -> [(assembly): list.first()] }

        config.species_name = referenceDict.find { k, v -> v instanceof Map && v.containsKey(config.reference) }?.key
        config.organism = config.species_name.toLowerCase().replace(' ', '_').split('\\(')[0].trim()
        config.organism_code = keggDict[config.species_name]

        // Add more logic here to set other config values based on the loaded information
        // This is a simplified version and would need to be expanded based on your specific requirements

        return config
    }

    static def loadMirna(config) {
        def referenceDict = jsonSlurper.parseText(new File("${config.globalResources}/reference_info/reference2.json").text)
        def keggDict = jsonSlurper.parseText(new File("${config.globalResources}/reference_info/kegg_reference.json").text)
        def organismTab = new File("${config.globalResources}/reference_info/organism_tab.tsv").readLines()
            .drop(1) // skip header
            .collect { line -> 
                def (assembly, full_name, kegg_term, release) = line.split('\t')
                [assembly: assembly, full_name: full_name, kegg_term: kegg_term, release: release]
            }
            .groupBy { it.assembly }
            .collectEntries { assembly, list -> [(assembly): list.first()] }

        config.species_name = referenceDict.find { k, v -> v instanceof Map && v.containsKey(config.reference) }?.key
        config.organism = config.species_name.toLowerCase().replace(' ', '_').split('\\(')[0].trim()
        config.organism_code = keggDict[config.species_name]

        // Add more logic here to set other config values based on the loaded information
        // This is a simplified version and would need to be expanded based on your specific requirements

        return config
    }

    static def loadAndConfigureUMI(config, wfConfigPath) {
        def wfConfig = jsonSlurper.parseText(new File(wfConfigPath).text)
        def primaryGuiParams = wfConfig.gui_params.primary
        def umiSettings = [:]
        def paramKeys = ["UMI_write_to", "UMI_R1_start", "UMI_R1_end", "insert_R1_start",
                         "UMI_R2_start", "UMI_R2_end", "insert_R2_start"]

        primaryGuiParams.UMI.list.each { umiType, _ ->
            def settings = [:]
            paramKeys.each { paramKey ->
                if (primaryGuiParams[paramKey]?.conditions?.value?.UMI?.containsKey(umiType)) {
                    settings[paramKey] = primaryGuiParams[paramKey].conditions.value.UMI[umiType]
                }
            }
            if (settings) {
                umiSettings[umiType] = settings
            }
        }

        def umiType = config.UMI
        if (umiSettings.containsKey(umiType)) {
            config.putAll(umiSettings[umiType])
        }

        return config
    }

    // Additional utility methods can be added here
}