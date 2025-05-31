"use client";

import type React from "react";

import { toPng } from "html-to-image";
import { useState, useRef } from "react";
import axios from "axios";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select";
import { Badge } from "@/components/ui/badge";
import { Alert, AlertDescription } from "@/components/ui/alert";
import {
  Upload,
  Download,
  Share,
  FileJson,
  FileSpreadsheet,
  AlertCircle,
  Dna,
} from "lucide-react";

export default function Home() {
  const [formData, setFormData] = useState({
    email: "",
    reference: "NM_007294.4",
    alignment_type: "global",
    variant_source: "clinvar",
    population: "",
  });
  const [file, setFile] = useState<File | null>(null);
  const [result, setResult] = useState<any>(null);
  const [loading, setLoading] = useState(false);

  const clinicalLabels = [
    "Pathogenic",
    "Likely pathogenic",
    "Benign",
    "Likely benign",
    "Uncertain significance",
    "Unknown",
  ];

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    if (!file) return;
    setLoading(true);
    const data = new FormData();
    Object.entries(formData).forEach(([key, value]) => data.append(key, value));
    data.append("file", file);

    try {
      const response = await axios.post("http://127.0.0.1:8000/analyze", data);
      setResult(response.data);
    } catch (err) {
      alert("Failed to analyze sequence.");
    } finally {
      setLoading(false);
    }
  };

  const handleFileUpload = (e: React.ChangeEvent<HTMLInputElement>) => {
    const selectedFile = e.target.files?.[0] || null;
    setFile(selectedFile);
  };

  const handleDragOver = (e: React.DragEvent) => {
    e.preventDefault();
  };

  const handleDrop = (e: React.DragEvent) => {
    e.preventDefault();
    const droppedFile = e.dataTransfer.files[0];
    if (droppedFile) {
      setFile(droppedFile);
    }
  };

  const exportJSON = () => {
    if (!result) return;
    const blob = new Blob([JSON.stringify(result, null, 2)], {
      type: "application/json",
    });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = `${result.sample_id}_report.json`;
    link.click();
    URL.revokeObjectURL(url);
  };

  const exportCSV = () => {
    if (!result || !result.variants) return;
    const header = Object.keys(result.variants[0]).join(",");
    const rows = result.variants
      .map((v: any) => Object.values(v).join(","))
      .join("\n");
    const csv = `${header}\n${rows}`;
    const blob = new Blob([csv], { type: "text/csv" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = `${result.sample_id}_report.csv`;
    link.click();
    URL.revokeObjectURL(url);
  };

  const downloadImage = async () => {
    if (!resultRef.current) return;
    const dataUrl = await toPng(resultRef.current);
    const link = document.createElement("a");
    link.href = dataUrl;
    link.download = `${result.sample_id}_report.png`;
    link.click();
  };

  const countVariantTypes = (type: string) => {
    return result?.variants?.filter((v: any) => v.type === type).length || 0;
  };

  const countClinicalSignificance = (label: string) => {
    return (
      result?.variants?.filter((v: any) => v.clinical_significance === label)
        .length || 0
    );
  };

  const getVariantTypeColor = (type: string) => {
    switch (type.toLowerCase()) {
      case "snp":
        return "bg-green-100 text-green-800 border-green-200";
      case "indel":
        return "bg-blue-100 text-blue-800 border-blue-200";
      default:
        return "bg-gray-100 text-gray-800 border-gray-200";
    }
  };

  const getClinicalSignificanceColor = (significance: string) => {
    switch (significance) {
      case "Pathogenic":
        return "bg-red-500 text-white";
      case "Likely pathogenic":
        return "bg-orange-500 text-white";
      case "Benign":
        return "bg-green-500 text-white";
      case "Likely benign":
        return "bg-green-400 text-white";
      case "Uncertain significance":
        return "bg-yellow-400 text-black";
      default:
        return "bg-gray-400 text-white";
    }
  };

  const getPathogenicVariants = () => {
    if (!result?.variants) return [];
    return result.variants.filter(
      (v: any) =>
        v.clinical_significance === "Pathogenic" ||
        v.clinical_significance === "Likely pathogenic"
    );
  };

  const resultRef = useRef<HTMLDivElement>(null);

  return (
    <div className="min-h-screen bg-gray-50 p-6">
      {/* Header */}
      <div className="max-w-7xl mx-auto mb-8">
        <div className="flex items-center justify-center gap-3 mb-2">
          <div className="w-8 h-8 bg-emerald-500 rounded-lg flex items-center justify-center">
            <Dna className="w-5 h-5 text-white" />
          </div>
          <h1 className="text-3xl font-bold text-gray-900">
            BRCA1 Sequence Analyzer
          </h1>
        </div>
      </div>

      <div className="max-w-7xl mx-auto grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Left Panel - Analysis Parameters */}
        <Card className="h-fit">
          <CardHeader>
            <CardTitle className="text-xl font-semibold text-gray-900">
              Analysis Parameters
            </CardTitle>
            <p className="text-sm text-gray-600">
              Configure the parameters for your BRCA1 sequence analysis
            </p>
          </CardHeader>
          <CardContent className="space-y-6">
            <form onSubmit={handleSubmit} className="space-y-3">
              <div className="space-y-2">
                <Label
                  htmlFor="email"
                  className="text-sm font-medium text-gray-700"
                >
                  Email (for NCBI Entrez queries)
                </Label>
                <Input
                  id="email"
                  type="email"
                  placeholder="your.email@example.com"
                  value={formData.email}
                  onChange={(e) =>
                    setFormData({ ...formData, email: e.target.value })
                  }
                  required
                  className="w-full"
                />
              </div>

              <div className="space-y-2">
                <Label
                  htmlFor="reference"
                  className="text-sm font-medium text-gray-700"
                >
                  Reference Sequence Accession
                </Label>
                <Input
                  id="reference"
                  type="text"
                  value={formData.reference}
                  onChange={(e) =>
                    setFormData({ ...formData, reference: e.target.value })
                  }
                  className="w-full"
                />
              </div>

              <div className="grid grid-cols-2 gap-4">
                <div className="space-y-2">
                  <Label className="text-sm font-medium text-gray-700">
                    Alignment Type
                  </Label>
                  <Select
                    value={formData.alignment_type}
                    onValueChange={(value) =>
                      setFormData({ ...formData, alignment_type: value })
                    }
                  >
                    <SelectTrigger>
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="global">Global</SelectItem>
                      <SelectItem value="local">Local</SelectItem>
                    </SelectContent>
                  </Select>
                </div>

                <div className="space-y-2">
                  <Label className="text-sm font-medium text-gray-700">
                    Variant Annotation Source
                  </Label>
                  <Select
                    value={formData.variant_source}
                    onValueChange={(value) =>
                      setFormData({ ...formData, variant_source: value })
                    }
                  >
                    <SelectTrigger>
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="clinvar">ClinVar</SelectItem>
                      <SelectItem value="none">None</SelectItem>
                    </SelectContent>
                  </Select>
                </div>
              </div>

              <div className="space-y-2">
                <Label
                  htmlFor="population"
                  className="text-sm font-medium text-gray-700"
                >
                  Population (optional)
                </Label>
                <Input
                  id="population"
                  type="text"
                  placeholder="e.g., European, African, Asian"
                  value={formData.population}
                  onChange={(e) =>
                    setFormData({ ...formData, population: e.target.value })
                  }
                  className="w-full"
                />
              </div>

              {/* File Upload Area */}
              <div className="space-y-2">
                <div className="flex justify-between items-center">
                  <Label className="text-sm font-medium text-gray-700">
                    Upload FASTA File
                  </Label>
                  <Button variant="outline" size="sm" type="button">
                    Paste Sequence
                  </Button>
                </div>
                <div
                  className="border-2 border-dashed border-gray-300 rounded-lg p-4 text-center hover:border-gray-400 transition-colors"
                  onDragOver={handleDragOver}
                  onDrop={handleDrop}
                >
                  <Upload className="w-6 h-6 text-gray-400 mx-auto mb-4" />
                  <p className="text-sm text-gray-600 mb-2">
                    Drag and drop your FASTA file here, or
                  </p>
                  <Button
                    variant="outline"
                    size="sm"
                    type="button"
                    onClick={() =>
                      document.getElementById("file-input")?.click()
                    }
                  >
                    Browse Files
                  </Button>
                  <input
                    id="file-input"
                    type="file"
                    className="hidden"
                    onChange={handleFileUpload}
                    accept=".fasta,.fa,.fas"
                  />
                  {file && (
                    <p className="text-sm text-green-600 mt-2">
                      Selected: {file.name}
                    </p>
                  )}
                </div>
              </div>

              <Button
                type="submit"
                className="w-full bg-emerald-500 hover:bg-emerald-600 text-white"
                disabled={loading || !file}
              >
                {loading ? "Analyzing..." : "Analyze Sequence"}
              </Button>
            </form>
          </CardContent>
        </Card>

        {/* Right Panel - Analysis Results */}
        {result && (
          <Card ref={resultRef}>
            <CardHeader>
              <div className="flex justify-between items-start">
                <div>
                  <CardTitle className="text-xl font-semibold text-gray-900">
                    Analysis Results
                  </CardTitle>
                  <p className="text-sm text-gray-600">
                    Detailed analysis of sample {result.sample_id}
                  </p>
                </div>
                <div className="flex gap-2">
                  <Button variant="outline" size="sm" onClick={exportJSON}>
                    <FileJson className="w-4 h-4 mr-1" /> JSON
                  </Button>
                  <Button variant="outline" size="sm" onClick={exportCSV}>
                    <FileSpreadsheet className="w-4 h-4 mr-1" />
                    CSV
                  </Button>
                  <Button variant="outline" size="sm" onClick={downloadImage}>
                    <Share className="w-4 h-4 mr-1" />
                    Share
                  </Button>
                  <Button
                    size="sm"
                    className="bg-emerald-500 hover:bg-emerald-600"
                    onClick={downloadImage}
                  >
                    <Download className="w-4 h-4 mr-1" />
                    Download Report
                  </Button>
                </div>
              </div>
            </CardHeader>
            <CardContent className="space-y-6">
              {/* Sample Information Grid */}
              <div className="grid grid-cols-3 gap-4 text-sm">
                <div>
                  <p className="text-gray-600">Sample ID</p>
                  <p className="font-semibold">{result.sample_id}</p>
                </div>
                <div>
                  <p className="text-gray-600">Reference ID</p>
                  <p className="font-semibold">{result.reference_id}</p>
                </div>
                <div>
                  <p className="text-gray-600">Sample Length</p>
                  <p className="font-semibold">{result.sample_length} bp</p>
                </div>
                <div>
                  <p className="text-gray-600">Alignment Score</p>
                  <p className="font-semibold">{result.alignment_score}</p>
                </div>
                <div>
                  <p className="text-gray-600">Processing Time</p>
                  <p className="font-semibold">4.89 seconds</p>
                </div>
                <div>
                  <p className="text-gray-600">Analysis Date</p>
                  <p className="font-semibold">Jan 15, 2025</p>
                  <p className="text-xs text-gray-500">09:45 AM</p>
                </div>
              </div>

              {/* Variant Summary */}
              <div>
                <h3 className="font-semibold text-gray-900 mb-4">
                  Variant Summary
                </h3>
                <div className="grid grid-cols-2 gap-6">
                  <div>
                    <p className="text-sm text-gray-600 mb-3">Variant Types</p>
                    <div className="space-y-2">
                      <div className="flex items-center justify-between">
                        <div className="flex items-center gap-2">
                          <div className="w-2 h-2 bg-green-500 rounded-full"></div>
                          <span className="text-sm">SNP</span>
                        </div>
                        <span className="text-sm font-semibold">
                          {countVariantTypes("SNP")}
                        </span>
                      </div>
                      <div className="flex items-center justify-between">
                        <div className="flex items-center gap-2">
                          <div className="w-2 h-2 bg-blue-500 rounded-full"></div>
                          <span className="text-sm">Indel</span>
                        </div>
                        <span className="text-sm font-semibold">
                          {countVariantTypes("indel")}
                        </span>
                      </div>
                    </div>
                  </div>
                  <p className="text-sm text-gray-600 mb-3">
                    Clinical Significance
                  </p>
                  <div className="space-y-2">
                    {clinicalLabels.map((label) => {
                      const count = countClinicalSignificance(label);
                      if (count === 0) return null;

                      return (
                        <div
                          key={label}
                          className="flex items-center justify-between"
                        >
                          <Badge
                            className={`${getClinicalSignificanceColor(
                              label
                            )} text-xs`}
                          >
                            {label}
                          </Badge>
                          <span className="text-sm font-semibold">{count}</span>
                        </div>
                      );
                    })}
                  </div>
                </div>
              </div>

              {/* Key Findings */}
              {getPathogenicVariants().length > 0 && (
                <div>
                  <h3 className="font-semibold text-gray-900 mb-4">
                    Key Findings
                  </h3>
                  <Alert className="border-red-200 bg-red-50">
                    <AlertCircle className="h-4 w-4 text-red-600" />
                    <AlertDescription>
                      <div className="space-y-2">
                        <p className="font-semibold text-red-800">
                          Pathogenic Variants Detected
                        </p>
                        <p className="text-sm text-red-700">
                          This sample contains {getPathogenicVariants().length}{" "}
                          pathogenic or likely pathogenic variant(s) that may be
                          clinically significant.
                        </p>
                        <ul className="text-sm text-red-700 space-y-1">
                          {getPathogenicVariants().map(
                            (variant: any, idx: number) => (
                              <li key={idx}>
                                • {variant.notation} -{" "}
                                {variant.clinical_significance}
                              </li>
                            )
                          )}
                        </ul>
                      </div>
                    </AlertDescription>
                  </Alert>
                </div>
              )}
            </CardContent>
          </Card>
        )}
      </div>
    </div>
  );
}
