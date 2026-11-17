// =====================================================================
// FISH_Prominence_Hilfe.ijm
// Fiji / ImageJ 1.54t
//
// Hilfsmakro zum Festlegen der PROMINENCE (Find Maxima) fuer FISH-Punkte.



requires("1.54t");
setOption("BlackBackground", true);

var useComma = false;

// ---------- Hilfsfunktionen ----------
function fmt(x, dec) {
    s = d2s(x, dec);
    if (useComma) s = replace(s, "\\.", ",");
    return s;
}

function closeID(id) {
    selectImage(id);
    close();
}

function blankMask(w, h) {
    newImage("blank", "8-bit black", w, h, 1);
    return getImageID();
}

function makeRoiMask(roiIdx, w, h) {
    newImage("roiMask", "8-bit black", w, h, 1);
    id = getImageID();
    roiManager("select", roiIdx);
    setColor(255);
    fill();
    run("Select None");
    return id;
}

function orMask(aID, bID) {
    imageCalculator("Max create", aID, bID);
    newID = getImageID();
    closeID(aID);
    closeID(bID);
    return newID;
}

function andMask(aID, bID) {
    imageCalculator("AND create", aID, bID);
    return getImageID();
}

function subMask(aID, bID) {
    imageCalculator("Subtract create", aID, bID);
    return getImageID();
}

function countWhite(id) {
    selectImage(id);
    run("Select None");
    getHistogram(values, counts, 256);
    return counts[255];
}

// ---------- Einstellungen ----------
inDir   = getDirectory("Ordner mit den .nd2-Dateien waehlen");
outDir  = inDir + "Auswertung_FISH_Neu" + File.separator;
helpDir = outDir + "Prominence_Hilfe" + File.separator;
File.makeDirectory(outDir);
File.makeDirectory(helpDir);

Dialog.create("Prominence-Hilfe (FISH)");
Dialog.addString("Gen / Sonde (fuer Dateinamen):", "Esco2");
Dialog.addNumber("Kanal der Sonde (2 oder 3):", 2);
Dialog.addString("Prominence-Werte (kommagetrennt, max. 11):", "10,25,50,100,200,400");
Dialog.addNumber("Standard-Ausschnittsgroesse in Pixel (falls du kein Rechteck zeichnest):", 500);
Dialog.addNumber("Vergroesserung der Kacheln:", 2);
Dialog.addNumber("Nur jedes n-te Bild verwenden (1 = alle):", 1);
Dialog.addCheckbox("Punktzahl auch ueber die gesamte Auswertflaeche berechnen (dauert bei grossen Mosaiken laenger)", true);
Dialog.addCheckbox("Dezimaltrennzeichen Komma statt Punkt", false);
Dialog.show();

marker   = Dialog.getString();
promText = Dialog.getString();
chNo     = Dialog.getNumber();
cropSize = Dialog.getNumber();
mag      = Dialog.getNumber();
stepN    = Dialog.getNumber();
doTotal  = Dialog.getCheckbox();
useComma = Dialog.getCheckbox();
if (stepN < 1) stepN = 1;
if (mag < 1) mag = 1;

pvals = newArray(0);
parts = split(promText, ",");
for (i = 0; i < parts.length; i++) {
    v = parseFloat(String.trim(parts[i]));
    if (!isNaN(v) && pvals.length < 11) pvals = Array.concat(pvals, v);
}
if (pvals.length == 0) exit("Keine gueltigen Prominence-Werte eingegeben.");
nP = pvals.length;

allFiles = getFileList(inDir);
files = newArray(0);
for (i = 0; i < allFiles.length; i++) {
    if (endsWith(toLowerCase(allFiles[i]), ".nd2"))
        files = Array.concat(files, allFiles[i]);
}
if (files.length == 0) exit("Keine .nd2-Dateien im gewaehlten Ordner gefunden.");

outTab = "Bild\tProminence\tPunkte_im_Ausschnitt\tPunkte_gesamt_Auswertflaeche\tAuswertflaeche_um2\tPunkte_pro_1000um2\n";

run("Close All");
roiManager("reset");
run("ROI Manager...");

for (f = 0; f < files.length; f = f + stepN) {
    name = files[f];
    baseName = substring(name, 0, lastIndexOf(name, "."));
    print("(" + (f + 1) + "/" + files.length + ") " + name);

    run("Bio-Formats Importer", "open=[" + inDir + name + "] color_mode=Composite view=Hyperstack stack_order=XYCZT");
    origID = getImageID();
    getDimensions(w, h, nC, nZ, nT);
    getPixelSize(pxUnit, pixelW, pixelH, pixelD);

    if (nC < chNo) {
        print("   -> uebersprungen: Bild hat nur " + nC + " Kanaele.");
        run("Close All");
        continue;
    }
    if (nZ > 1) {
        run("Z Project...", "projection=[Max Intensity]");
        zID = getImageID();
        selectImage(origID);
        close();
        origID = zID;
    }
    if (nC > 1) Stack.setDisplayMode("composite");

    // ---------- Auswertbereich ----------
    roiPath = outDir + baseName + "_ROIs.zip";
    roiManager("reset");
    if (File.exists(roiPath)) {
        roiManager("open", roiPath);
        print("   -> vorhandene ROIs verwendet.");
    } else {
        selectImage(origID);
        run("Select None");
        setTool("freehand");
        waitForUser("Bereich einzeichnen: " + name,
            "Metastase (oder groben repraesentativen Bereich) umranden, dann OK.\n" +
            "Dieser Bereich wird NICHT gespeichert.");
        if (selectionType() != -1) {
            roiManager("add");
            roiManager("select", roiManager("count") - 1);
            roiManager("rename", "Meta_1");
            selectImage(origID);
            run("Select None");
        }
    }
    nTotal = roiManager("count");
    if (nTotal == 0) {
        print("   -> kein Bereich, uebersprungen.");
        run("Close All");
        continue;
    }

    // Auswertflaeche = Vereinigung der Metastasen minus Nekrose/Ausschluss
    metaU = blankMask(w, h);
    exclU = blankMask(w, h);
    for (i = 0; i < nTotal; i++) {
        selectImage(origID);
        roiManager("select", i);
        rname = Roi.getName();
        run("Select None");
        if (startsWith(rname, "Meta_")) {
            mk = makeRoiMask(i, w, h);
            metaU = orMask(metaU, mk);
        } else if (startsWith(rname, "Nekrose_") || startsWith(rname, "Ausschluss_")) {
            ek = makeRoiMask(i, w, h);
            exclU = orMask(exclU, ek);
        }
    }
    validID = subMask(metaU, exclU);
    closeID(metaU);
    closeID(exclU);
    validPx = countWhite(validID);
    if (validPx == 0) {
        print("   -> leere Auswertflaeche, uebersprungen.");
        run("Close All");
        continue;
    }
    validArea = validPx * pixelW * pixelH;

    // ---------- Kanal (ohne Kalibrierung, damit Koordinaten in Pixeln sind) ----------
    selectImage(origID);
    run("Select None");
    run("Duplicate...", "title=POI duplicate channels=" + chNo);
    chID = getImageID();
    run("Set Scale...", "distance=0 known=0 unit=pixel");

    // ---------- Ausschnitt waehlen ----------
    selectImage(origID);
    run("Select None");
    setTool("rectangle");
    waitForUser("Ausschnitt waehlen: " + name,
        "Zeichne ein RECHTECK um einen repraesentativen Ausschnitt mit FISH-Punkten\n" +
        "(am besten ein paar Dutzend Zellen breit, mit hellen und schwachen Punkten).\n" +
        "Kein Rechteck gezeichnet -> Standard-Ausschnitt in der Mitte der ersten Metastase.");
    cropW = cropSize;
    cropH = cropSize;
    if (selectionType() == 0) {
        getSelectionBounds(cx0, cy0, cropW, cropH);
    } else {
        selectImage(origID);
        roiManager("select", 0);
        getSelectionBounds(bx, by, bw, bh);
        run("Select None");
        cx0 = bx + bw / 2 - cropW / 2;
        cy0 = by + bh / 2 - cropH / 2;
    }
    selectImage(origID);
    run("Select None");
    cropW = minOf(cropW, w);
    cropH = minOf(cropH, h);
    cx0 = round(maxOf(0, minOf(cx0, w - cropW)));
    cy0 = round(maxOf(0, minOf(cy0, h - cropH)));

    selectImage(chID);
    run("Select None");
    makeRectangle(cx0, cy0, cropW, cropH);
    run("Duplicate...", "title=cropRaw");
    cropRawID = getImageID();
    run("Select None");

    // Anzeige-Grundbild: Kontrast gestreckt, RGB, vergroessert
    selectImage(cropRawID);
    run("Duplicate...", "title=cropDisp");
    run("Enhance Contrast", "saturated=0.35");
    run("8-bit");
    run("RGB Color");
    dispID = getImageID();
    run("Scale...", "x=" + mag + " y=" + mag + " interpolation=None create");
    baseID = getImageID();
    selectImage(dispID);
    close();

    nTiles = 1 + nP;
    tileIDs = newArray(nTiles);
    tileLabels = newArray(nTiles);
    totals = newArray(nP);

    selectImage(baseID);
    run("Select None");
    run("Duplicate...", "title=tile0");
    tileIDs[0] = getImageID();
    tileLabels[0] = "Rohbild (Ausschnitt)";

    rr = maxOf(4, round(mag * 3));

    for (t = 0; t < nP; t++) {
        pv = pvals[t];

        // Punkte im Ausschnitt
        selectImage(cropRawID);
        run("Select None");
        run("Duplicate...", "title=cropP");
        cpID = getImageID();
        run("Clear Results");
        run("Find Maxima...", "prominence=" + pv + " output=List");
        nS = nResults;
        xs = newArray(0);
        ys = newArray(0);
        if (nS > 0) {
            xs = Table.getColumn("X");
            ys = Table.getColumn("Y");
        }
        run("Clear Results");
        selectImage(cpID);
        close();

        selectImage(baseID);
        run("Select None");
        run("Duplicate...", "title=tileP");
        tID = getImageID();
        setColor("red");
        setLineWidth(2);
        for (s = 0; s < nS; s++) {
            drawOval(round(xs[s] * mag) - rr, round(ys[s] * mag) - rr, 2 * rr, 2 * rr);
        }
        tileIDs[1 + t] = tID;

        // Punkte ueber die gesamte Auswertflaeche
        totalCount = -1;
        perK = -1;
        if (doTotal) {
            selectImage(validID);
            run("Create Selection");
            selectImage(chID);
            run("Restore Selection");
            run("Clear Results");
            run("Find Maxima...", "prominence=" + pv + " output=Count");
            totalCount = getResult("Count", nResults - 1);
            run("Select None");
            run("Clear Results");
            perK = totalCount / (validArea / 1000);
            totals[t] = totalCount;
        }

        lbl = "Prominence " + pv + ": " + nS + " Punkte im Ausschnitt";
        if (doTotal) lbl = lbl + " | " + totalCount + " gesamt";
        tileLabels[1 + t] = lbl;

        totText = "-";
        perText = "-";
        if (doTotal) {
            totText = "" + totalCount;
            perText = fmt(perK, 3);
        }
        outTab = outTab + name + "\t" + fmt(pv, 1) + "\t" + nS + "\t" + totText + "\t" + fmt(validArea, 1) + "\t" + perText + "\n";
    }

    // ---------- Montage ----------
    selectImage(tileIDs[0]);
    tileW = getWidth();
    tileH = getHeight();
    cols = 3;
    rows = floor((nTiles + cols - 1) / cols);
    newImage("montage", "RGB black", cols * tileW, rows * tileH, 1);
    montID = getImageID();
    fs = 24;
    for (i = 0; i < nTiles; i++) {
        px = (i % cols) * tileW;
        py = floor(i / cols) * tileH;
        selectImage(tileIDs[i]);
        run("Select All");
        run("Copy");
        selectImage(montID);
        makeRectangle(px, py, tileW, tileH);
        run("Paste");
        run("Select None");
        setFont("SansSerif", fs, "bold");
        tw = getStringWidth(tileLabels[i]);
        setColor("black");
        fillRect(px + 4, py + 4, tw + 16, fs + 14);
        setColor("white");
        drawString(tileLabels[i], px + 12, py + fs + 8);
    }
    legend = "Montage " + baseName + " (" + marker + ")\nAusschnitt: x=" + cx0 + ", y=" + cy0 + ", " + cropW + " x " + cropH + " Pixel\n";
    for (i = 0; i < nTiles; i++) {
        legend = legend + "Zeile " + (floor(i / cols) + 1) + ", Spalte " + ((i % cols) + 1) + ": " + tileLabels[i] + "\n";
    }
    File.saveString(legend, helpDir + "Legende_" + baseName + "_" + marker + ".txt");
    selectImage(montID);
    saveAs("PNG", helpDir + "Montage_" + baseName + "_" + marker + ".png");

    // ---------- Kurve Punktzahl vs. Prominence ----------
    if (doTotal && nP > 1) {
        Plot.create("Punktzahl vs. Prominence: " + baseName, "Prominence", "Punkte gesamt (Auswertflaeche)");
        Plot.setLogScaleX();
        Plot.add("line", pvals, totals);
        Plot.add("circles", pvals, totals);
        Plot.show();
        saveAs("PNG", helpDir + "Kurve_" + baseName + "_" + marker + ".png");
    }

    File.saveString(outTab, helpDir + "Prominence_" + marker + ".txt");

    run("Clear Results");
    run("Close All");
    roiManager("reset");
}

print("Fertig. Montagen, Kurven und Tabelle in: " + helpDir);
showMessage("Fertig", "Prominence-Hilfe abgeschlossen.\n\nMontagen, Kurven und Tabelle liegen in:\n" + helpDir);
